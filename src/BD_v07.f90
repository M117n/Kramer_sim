!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%***************************************************************************%%
!%%**  PROGRAM       BROWNIAN DYNAMICS                                      **%%
!%%**  AUTHOR        ALPIXELS, MARTIN MNARIQUE                              **%%
!%%**  LICENSE       LGPL-V3                                                **%%
!%%**                                                                       **%%
!%%**  ENSEMBLE      NVT-VELOCITY SCALING                                   **%%
!%%**  ALGORITHM     ERMACK & MACCAMON                                      **%%
!%%**  DATE          julio 16, 2015                                         **%%
!%%**  VERSION       0.0, BETA TEST                                         **%%
!%%**                                                                       **%%
!%%**  OBSERVATIONS                                                         **%%
!%%**                                                                       **%%
!%%***************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%   MODULE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE BDVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: D       = KIND(1.0D0) !PROGRAM PRECISION
 INTEGER, PARAMETER:: MP      = 10000        !MAXIMUM NUMBER OF PARTICLES
 INTEGER, PARAMETER:: OUTDISP = 500
 INTEGER, PARAMETER:: OUTFILM = 500
 INTEGER, PARAMETER:: BDSTEP  = 500000       !SIMULATION STEPS
 INTEGER, PARAMETER:: NHILOS  = 2            !NUMERO DE HILOS A USAR
 INTEGER, PARAMETER:: CHUNK   = 50           !DYNAMIC ASSIGNATION
 INTEGER, PARAMETER:: UPDATELIST = 10
 INTEGER, PARAMETER:: NPAIRX = 1000000      !MAX NUMBER NEIGHBOR PAIRS
 INTEGER, PARAMETER:: NC = 4
 INTEGER, PARAMETER:: NPART=4*NC**3         !NUMBER OF PARTICLES
 
!COORRELACION (MSD)
 INTEGER, PARAMETER:: OUTNSAM = 10          !CORRELATION TIME
 INTEGER, PARAMETER:: TMAXD   = 100000      !MAXIMUM NUMBER OF INSTANT TIMES
 INTEGER, PARAMETER:: IT0     = 10          !TAKE A NEW t=0
 INTEGER, PARAMETER:: T0MAX   = 200000      !MAXIMUM NUMBER OF t=0 
 INTEGER, PARAMETER:: NSAMPLE = T0MAX*IT0    !NUMERO MAXIMO DE MUESTRAS A TOMAR
 INTEGER:: NTIME(NSAMPLE), TIME0(NSAMPLE)
 INTEGER:: TMAX,NTEL,T0
 REAL(D), PARAMETER:: TOBS   = 30.0         !OBSERVATION TIME
 REAL(D):: RX0(NPART,TMAXD),RY0(NPART,TMAXD),RZ0(NPART,TMAXD)
 REAL(D):: MSDX(TMAXD,NPART),MSDY(TMAXD,NPART),MSDZ(TMAXD,NPART),MSDT(TMAXD,NPART)
 REAL(D):: DTIME
 !****
 
 !GdR
 INTEGER, PARAMETER:: NHIS=???
 INTEGER:: SWITCH1
 REAL(D):: G(NHIS)



REAL(D), PARAMETER:: PI      = 4.D0*DATAN(1.D0) 
 REAL(D), PARAMETER:: SIGMA12 = 1.0D0             !DIAMETRO DE LA PARTICULA
 REAL(D), PARAMETER:: ESTAR12 = 1.0D0

 INTEGER:: IDUMM,IBD,NPAIR,NPARA
 INTEGER:: IFRAME 
 INTEGER:: SWITCH
 INTEGER:: VLI(NPAIRX), VLJ(NPAIRX)

 REAL(D):: RX(MP),RY(MP),RZ(MP)
 REAL(D):: FX(MP),FY(MP),FZ(MP)
 REAL(D):: RHOSTAR,PHI,EPOT,RCUT,RV
 REAL(D):: BOXX,BOXY,BOXZ,IBOXX,IBOXY,IBOXZ
 REAL(D):: A2,DK,H,DTT


END MODULE BDVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM BROWNIAN_DYNAMICS
 USE OMP_LIB
 USE BDVAR
 IMPLICIT NONE
 !$CALL OMP_SET_NUM_THREADS(NHILOS)

 OPEN(UNIT=13,FILE='MovieBD.xyz')

 CALL SETUP                             !DEFINE SIMULATION PARAMETERS
 CALL CONFCC                            !SIMULATION INITIAL CONFIGURATION

 CALL VERLET                            !VERLET NEIGHBORHOOD LIST
 CALL BDFORCE                           !COMPUTE FORCE
 CALL SINGI
 SWITCH=0                               !INICIALIZACION PARA CALCULAR FUNCION DE COORRELACION (MSD)
 SWITCH1=0                              !INICIALIZACION PARA CALCULAR GdR
 CALL DIFUSION 
 CALL GdR 

 DO IBD=1,BDSTEP
    CALL INTEGRATE  
    CALL BDFORCE
    CALL DYNAMICS
    IF (MOD(IBD,UPDATELIST) .EQ. 0) CALL VERLET    !ACTUALIZAR LA LISTA DE VERLET CADA UPDATELIST-NUMERO DE PASOS
    IF (MOD(IBD,OUTNSAM) .EQ. 0) CALL DIFUSION     !COMENZAR TOMANDO MUESTRAS PARA MSD
    IF (MOD(IBD,GDRPAST) .EQ. 0) CALL GdR          !TOMAR MUESTRAS PARA CALCULAR GdR
 ENDDO

  SWITCH=2           !CREAR LA ESTADISTICA DE MSD
  SWITCH1=2          !ESTADISTICA PARA GdR
  CALL DIFUSION

 CLOSE(UNIT=13)

 STOP
END PROGRAM BROWNIAN_DYNAMICS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%   SUBROUTINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE SETUP
 USE BDVAR
 IMPLICIT NONE
 REAL(D):: LBOX

 PHI=0.24                                   !PACKING FRACTION
 RHOSTAR=6.0D0*PHI/PI                       !SYSTEM DENSITY
 BOXX=(DBLE(NPART)/RHOSTAR)**(1.D0/3.D0)    !SIMULATION BOX SIDE LENGHT
 BOXY=BOXX
 BOXZ=BOXX
 IBOXX=1.D0/BOXX                            !INVERSE SIMULATION BOX SIDE LENGHT
 IBOXY=1.D0/BOXY
 IBOXZ=1.D0/BOXZ
 IDUMM=123456789                             !SEED FOR RANDON NUMBER GENERATOR
 LBOX=DMAX1(BOXX,BOXY,BOXZ)
 RCUT=LBOX/2.0D0                            !CUT-OFF 
 DK=1.5281
 A2=2450*(DEXP(DK/2.0D0)/(1.0D0 + 0.5D0*DK))**2
 DTT=0.0001                                 !STEP TIME TO THERMALIZE
 H=DTT
 IFRAME=0

 !DIFUSION
 T0=0
 NTEL=0
 DTIME=H*DBLE(OUTNSAM)                    !CORRELATION TIME STEP
 TMAX=INT(TOBS/DTIME)                     !TOTAL NUMBER OF TIME STEP

 RETURN
END SUBROUTINE SETUP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE CONFIG
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I,J
 REAL(D):: RAN1,X,Y,Z,DX,DY,DZ,RIJ
        
 X=RAN1(IDUMM)
 RX(1)=(X-0.5)*BOXX
 Y=RAN1(IDUMM)
 RY(1)=(Y-0.5)*BOXY
 Z=RAN1(IDUMM)
 RZ(1)=(Z-0.5)*BOXZ

 DO I=2,NPART
   1 CONTINUE
   X=RAN1(IDUMM)
   X=(X-0.5)*BOXX
   Y=RAN1(IDUMM)
   Y=(Y-0.5)*BOXY
   Z=RAN1(IDUMM)
   Z=(Z-0.5)*BOXZ

   DO J=1,I-1
      DX=X-RX(J)
      DY=Y-RY(J)
      DZ=Z-RZ(J)

      DX=DX-BOXX*ANINT(DX*IBOXX)
      DY=DY-BOXY*ANINT(DY*IBOXY)
      DZ=DZ-BOXZ*ANINT(DZ*IBOXZ)

      RIJ=DSQRT(DX*DX+DY*DY+DZ*DZ)
      IF(RIJ .LE. 1.D0) GOTO 1
   ENDDO

   RX(I)=X
   RY(I)=Y
   RZ(I)=Z
 ENDDO
  

 OPEN(UNIT=11,FILE='PicBD.xyz')
 WRITE(11,*)NPART
 WRITE(11,*)'FRAME',1

 DO I=1,NPART
    WRITE(11,*)'C',RX(I),RY(I),RZ(I)
 ENDDO

 RETURN
END SUBROUTINE CONFIG
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE CONFCC
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: IX,IY,IZ,LAT,M,I,K
 REAL(D):: CELLX,CELLY,CELLZ
 REAL(D):: CELLUX,CELLUY,CELLUZ
 REAL(D):: DX,DY,DZ,RIJ

 M=0
 CELLX=BOXX/DBLE(NC)                        !UNITARY CELL X LENGHT
 CELLY=BOXY/DBLE(NC)                        !UNITARY CELL Y LENGHT
 CELLZ=BOXZ/DBLE(NC)                        !UNITARY CELL Z LENGHT

 CELLUX=CELLX*0.5D0                         !CENTER OF CELL X DIRECTION
 CELLUY=CELLY*0.5D0                         !CENTER OF CELL Y DIRECTION
 CELLUZ=CELLZ*0.5D0                         !CENTER OF CELL Z DIRECTION

 RX(1)=0.0D0                                !X POSITION ONE
 RY(1)=0.0D0                                !Y POSITION ONE
 RZ(1)=0.0D0                                !Z POSITION ONE

 RX(2)=CELLUX                               !X POSITION TWO
 RY(2)=CELLUY                               !Y POSITION TWO
 RZ(2)=0.0D0                                !Z POSITION TWO

 RX(3)=0.0D0                                !X POSITION THREE
 RY(3)=CELLUY                               !Y POSITION THREE
 RZ(3)=CELLUZ                               !Z POSITION THREE

 RX(4)=CELLUX                               !X POSITION FOUR
 RY(4)=0.0D0                                !Y POSITION FOUR
 RZ(4)=CELLUZ                               !Z POSITION FOUR


 DO IX=1,NC
    DO IY=1,NC
       DO IZ=1,NC
          DO LAT=1,4
             RX(LAT+M)=RX(LAT)+CELLX*DBLE(IX-1)
             RY(LAT+M)=RY(LAT)+CELLY*DBLE(IY-1)
             RZ(LAT+M)=RZ(LAT)+CELLZ*DBLE(IZ-1)
          ENDDO
          M=M+4
       ENDDO
    ENDDO
 ENDDO

 DO I=1,NPART                               !POSITIONS AT BOX CENTER
    RX(I)=RX(I)-0.5D0*BOXX
    RY(I)=RY(I)-0.5D0*BOXY
    RZ(I)=RZ(I)-0.5D0*BOXZ
 ENDDO

 RETURN
END SUBROUTINE CONFCC
!*******************************************************************************************************************
!*******************************************************************************************************************
!*******************************************************************************************************************
 SUBROUTINE VERLET
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I,J
 REAL(D):: DX,DY,DZ,RIJ

 NPAIR=0

 RV= RCUT + 0.3
 DO I=1, NPART -1 
   DO J=I+1, NPART
     DX=RX(I)-RX(J)
     DY=RY(I)-RY(J)
     DZ=RZ(I)-RZ(J)

     DX=DX-BOXX*ANINT(DX*IBOXX)
     DY=DY-BOXY*ANINT(DY*IBOXY)
     DZ=DZ-BOXZ*ANINT(DZ*IBOXZ)
		
     RIJ=DSQRT(DX*DX+DY*DY+DZ*DZ)
     IF(RIJ .LE. RV) THEN
        NPAIR=NPAIR + 1
        VLI(NPAIR)=I
        VLJ(NPAIR)=J
     END IF
   END DO
 END DO


END SUBROUTINE VERLET
!********************************************************************************
!********************************************************************************
!********************************************************************************
SUBROUTINE SINGI
 USE BDVAR
 IMPLICIT NONE

 OPEN(UNIT=12,FILE='WazapBD.txt')
 OPEN(UNIT=77,FILE='EPOT.txt')                 !IMPRIMIR NUMERO DE PASOS CON ENERGIA POTENCIAL
 !OPEN(UNIT=78,FILE='POSICIONES.txt')           !IMPRIMIR POSICIONES

 EPOT=EPOT/DBLE(NPART)

 WRITE(6,01) RHOSTAR,NPART,PHI,BOXX,EPOT,BDSTEP
 WRITE(12,01) RHOSTAR,NPART,PHI,BOXX,EPOT,BDSTEP

 CLOSE(UNIT=12)

 01 FORMAT(1X,//,'***  SIMULATION PARAMETERS  *** ',// &
          'REDUCED DENSITY               ',F15.8,/ &
          'NUMBER OF PARTICLES           ',I10  ,/ &
          'PACKING                       ',F15.8,/ &
          'SIMULATION X BOX SIZE         ',F15.8,/ &
          'POTENTIAL ENERGY              ',F15.8,/ &
          'SIMULATION STEPS              ',I10  ,//)

 RETURN
END SUBROUTINE SINGI
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE BDFORCE
 USE BDVAR
 IMPLICIT NONE
 REAL(D), PARAMETER:: RCW = 1.16499305075
 INTEGER:: I,J,K
 REAL(D):: X,Y,Z,DX,DY,DZ,RIJ,RIJSQ
 REAL(D):: BUR,BETAUR,FZA,FORCE
 REAL(D):: BURWL,BETAURWL,FORCEWL
 REAL(D):: ECOL,EWAL,DZW
 
 !falta la interacción de la particula np con las paredes, modificar

 EPOT=0.0D0
 EWAL=0.0D0 
 ECOL=0.0D0

 DO I=1,NPART
    FX(I)=0.0D0
    FY(I)=0.0D0
    FZ(I)=0.0D0
 END DO
 DO K=1,NPAIR      !SE CONTARA LA FUERZA POR PARES DE VECINOS
    I=VLI(K)
    J=VLJ(K) 

    DX=RX(I)-RX(J)
    DY=RY(I)-RY(J)
    DZ=RZ(I)-RZ(J)

    DX=DX - ANINT(DX*IBOXX)*BOXX
    DY=DY - ANINT(DY*IBOXY)*BOXY
    DZ=DZ - ANINT(DZ*IBOXZ)*BOXZ

    RIJSQ=DX*DX + DY*DY + DZ*DZ
    RIJ=DSQRT(RIJSQ)

    IF (RIJ .LT. RCUT) THEN
       BUR=BETAUR(RIJ,A2,DK)
       ECOL=ECOL + BUR
       !PRINT*, ECOL
         
       FZA=FORCE(RIJ,RIJSQ,A2,DK)
       FX(I)=FX(I) + FZA*DX
       FY(I)=FY(I) + FZA*DY
       FZ(I)=FZ(I) + FZA*DZ

       FX(J)=FX(J) - FZA*DX
       FY(J)=FY(J) - FZA*DY
       FZ(J)=FZ(J) - FZA*DZ
     ENDIF
 END DO 

 EPOT=ECOL + EWAL

 RETURN
 END SUBROUTINE BDFORCE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INTEGRATE
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(D):: SIGMA,RDX,RDY,RDZ,GASDEV
 REAL(D):: RXN,RYN,RZN,TEMP

 SIGMA=DSQRT(2.0D0*H)
 TEMP=0
 DO I=1,NPART
    RDX=SIGMA*GASDEV(IDUMM)
    RDY=SIGMA*GASDEV(IDUMM)
    RDZ=SIGMA*GASDEV(IDUMM)
    TEMP=TEMP + H
    !WRITE(78,*) RX(I),RY(I),RZ(I), TEMP

    RXN=RX(I) + FX(I)*H + RDX
    RYN=RY(I) + FY(I)*H + RDY
    RZN=RZ(I) + FZ(I)*H + RDZ

    RXN=RXN - ANINT(RXN*IBOXX)*BOXX
    RYN=RYN - ANINT(RYN*IBOXY)*BOXY
    RZN=RZN - ANINT(RZN*IBOXZ)*BOXZ

    RX(I)=RXN
    RY(I)=RYN
    RZ(I)=RZN
 ENDDO 

 RETURN
END SUBROUTINE INTEGRATE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE DYNAMICS
 USE BDVAR
 IMPLICIT NONE

 EPOT=EPOT/DBLE(NPART)

 !IF(MOD(IBD,OUTDISP) .EQ. 0)THEN        !MANDAR A IMPRIMIR A PANTALLA EL NUMERO DE PASOS CON LA ENERGIA POTENCIAL
   !WRITE(6,03)IBD,EPOT
 !ENDIF

IF(MOD(IBD,500) .EQ. 0)THEN
	WRITE(77,*) IBD,EPOT
END IF

 IF((MOD(IBD,OUTFILM) .EQ. 0) .AND. IBD .GE. 75000)CALL MOVIE

 SWITCH=1              !PARA CALCULAR LAS FUNCIONES DE COORRELACION
 SWITCH1=1             !PARA TOMAR MUESTRAS DE GdR

 03 FORMAT(I10,F10.5)

 RETURN
END SUBROUTINE DYNAMICS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE DIFUSION
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I,T,K,TT0,DELT
 REAL(D):: CTIME

 
 IF(SWITCH .EQ. 0)THEN                      !SET TO ZERO CORRELATION FUNCTIONS
   DO I=1,TMAX
      NTIME(I)=0                            !NUMBER OF SAMPLE FOR TIME ti
      DO K=1,npart                              !COLLOID & SOLVENT <A(0)A(ðt)>
         MSDX(I,K)=0.0D0                      !X MEAN SQUARE DISPLACEMENT
         MSDY(I,K)=0.0D0                      !Y MEAN SQUARE DISPLACEMENT
         MSDZ(I,K)=0.0D0                      !Z MEAN SQUARE DISPLACEMENT
         MSDT(I,K)=0.0D0                      !TOTAL MEAN SQUARE DISPLACEMENT
      ENDDO
   ENDDO

 ELSEIF(SWITCH .EQ. 1)THEN                  !COMPUTE CORRELATION FUNCTIONS
   NTEL=NTEL+1

   IF(MOD(NTEL,IT0) .EQ. 0)THEN             !TAKE A NEW t=0
     T0=T0+1                                !NUMBER OF t=0
     TT0=MOD(T0-1,T0MAX) + 1
     TIME0(TT0)=NTEL
     DO I=1,npart
        
        RX0(I,TT0)=RX(I)                  !RX(t=0)
        RY0(I,TT0)=RY(I)                  !RY(t=0)
        RZ0(I,TT0)=RZ(I)                  !RZ(t=0)

     ENDDO    
   ENDIF

   DO T=1,MIN(T0,T0MAX)
      DELT=NTEL - TIME0(T) + 1
      IF(DELT .LE. TMAX)THEN
        NTIME(DELT)=NTIME(DELT) + 1
        DO I=1,NPART                               !***********
           K=NPART+1 
           !<[x(0) - x(t)]^2>
           MSDX(DELT,I)=MSDX(DELT,I) + (RX(I) - RX0(I,T))**2     !*****************
           MSDY(DELT,I)=MSDY(DELT,I) + (RY(I) - RY0(I,T))**2
           MSDZ(DELT,I)=MSDZ(DELT,I) + (RZ(I) - RZ0(I,T))**2
           MSDT(DELT,I)=MSDT(DELT,I) + (RX(I) - RX0(I,T))**2 &
                                     + (RY(I) - RY0(I,T))**2 &
                                     + (RZ(I) - RZ0(I,T))**2     !**********************************
        ENDDO
      ENDIF
   ENDDO
 ELSEIF(SWITCH .EQ. 2)THEN                  !COMPUTE AVERAGES
   OPEN(UNIT=15,FILE='MSD.dat')
   
   DO I=1,TMAX                              !AVERAGE OF CORRELATION FUNCTIONS
      DO K=1,NPART
         MSDX(I,K)=MSDX(I,K)/DBLE(NTIME(I)*NPART)
         MSDY(I,K)=MSDY(I,K)/DBLE(NTIME(I)*NPART)
         MSDZ(I,K)=MSDZ(I,K)/DBLE(NTIME(I)*NPART)
         MSDT(I,K)=MSDT(I,K)/DBLE(NTIME(I)*NPART)
      ENDDO
      CTIME=H*REAL(OUTNSAM)*(REAL(I))
      WRITE(15,*)CTIME,MSDX(I,K),MSDY(I,K),MSDZ(I,K),MSDT(I,K)
   ENDDO

   CLOSE(UNIT=15)
 ENDIF

 RETURN
END SUBROUTINE DIFUSION
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
SUBROUTINE GdR

 IMPLICIT NONE
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: NGR, IG
 REAL(D):: DELG

 IF (SWITCH1.EQ.0) THEN      !INICIALIZACION
    NGR=0
    DELG= BOXX/(2*NHIS)
    DO I=1, NHIS
       G(I)=0
    END DO
 END IF

 IF (SWITCH1.EQ.1) THEN      !COMENZAR A TOMAR MUESTRAS
    NGR=NGR+1
    DO I=1, NPART -1
       DO J=1, NPART
          XR=X(I) - X(J) 
          YR=Y(I) - Y(J)
          ZR=Z(I) - Z(J)
  
          XR=XR - ANINT(RXN*IBOXX)*BOXX     !CONDICIONES DE FRONTERA
          YR=YR - ANINT(RYN*IBOXY)*BOXY
          ZR=ZR - ANINT(RZN*IBOXZ)*BOXZ
 
          R=DSQRT(XR*XR+YR*YR+ZR*ZR)
          IF (R.LT.(BOXX/2)) THEN             !RADIO A LA MITAD DE LA CAJA
             IG=INT(R/DELG)
             G(IG)=G(IG) + 2
          END IF
       END DO
    END DO
 END IF
 
 IF (SWITCH.EQ.2) THEN                   !CALCULAR GdR (ESTADISTICA)
    DO I=1, NHIS
       R=DELG*(I+0.5)
       VB=((I+1)**3-I**3)*DELG**3              !VOLUMEN DONDE CALCULAREMOS GdR
       NID=(4/3)*PI*VB*RHOSTAR
       G(I)=G(I)/(NGR*NPART*NID)
    END DO
 END IF




END SOBROUTINE GdR
!********************************************************************************
!********************************************************************************
!********************************************************************************
 SUBROUTINE MOVIE
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I

 IFRAME=IFRAME + 1
 
 WRITE(13,*)NPART
 WRITE(13,*)'FRAME',IFRAME

 DO I=1,NPART
    WRITE(13,*)'C',RX(I),RY(I),RZ(I)
 ENDDO

 RETURN
END SUBROUTINE MOVIE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION RAN1(IDUMM)
 INTEGER, PARAMETER:: D = KIND(1.0D0)
 INTEGER:: IDUMM,IA,IM,IQ,IR,NTAB,NDIV
 REAL(D):: RAN1,AM,EPS,RNMX
 PARAMETER( IA=16807, IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
            NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.0-EPS )
 INTEGER:: J,K,IV(NTAB),IY
 SAVE IV,IY
 DATA IV /NTAB*0/,IY /0/

 IF(IDUMM .LE. 0 .OR.  IY .EQ. 0)THEN
   IDUMM=MAX(-IDUMM,1)
   DO J=NTAB+8,1,-1
      K=IDUMM/IQ
      IDUMM=IA*(IDUMM-K*IQ)-IR*K
      IF(IDUMM .LT. 0)IDUMM=IDUMM + IM
      IF(J .LE. NTAB)IV(J)=IDUMM
   ENDDO
   IY=IV(1)
 ENDIF

 K=IDUMM/IQ
 IDUMM=IA*(IDUMM-K*IQ) - IR*K
 IF(IDUMM .LT. 0)IDUMM=IDUMM + IM
 J=1 + IY/NDIV
 IY=IV(J)
 IV(J)=IDUMM
 RAN1=MIN(AM*IY,RNMX)

 RETURN
END FUNCTION RAN1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION GASDEV(IDUMM)
 INTEGER, PARAMETER:: D = KIND(1.0D0)
 INTEGER:: IDUMM
 REAL(D):: GASDEV
 INTEGER:: ISET
 REAL(D):: FAC,GSET,RSQ,V1,V2,RAN1
 SAVE ISET,GSET
 DATA ISET/0/


 IF(IDUMM .LT. 0)ISET=0
 IF(ISET .EQ. 0)THEN
   13 CONTINUE
   V1=2.0D0*RAN1(IDUMM) - 1.0D0
   V2=2.0D0*RAN1(IDUMM) - 1.0D0
   RSQ=V1**2 + V2**2
   IF( RSQ .GE. 1.0D0 .OR. RSQ .EQ. 0.0D0)GOTO 13
   FAC=DSQRT(-2.0*DLOG(RSQ)/RSQ)
   GSET=V1*FAC
   GASDEV=V2*FAC
   ISET=1
 ELSE
   GASDEV=GSET
   ISET=0
 ENDIF

 RETURN
END FUNCTION GASDEV
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION BETAUR(RIJ,A2,DK)
 IMPLICIT NONE
 INTEGER, PARAMETER:: D  = KIND(1.0D0)
 REAL(D):: BETAUR,RIJ,A2,DK

 BETAUR=A2*DEXP(-DK*RIJ)/RIJ
 
 RETURN
END FUNCTION BETAUR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION FORCE(RIJ,RIJSQ,A2,DK)
 IMPLICIT NONE
 INTEGER, PARAMETER:: D  = KIND(1.0D0)
 REAL(D):: FORCE,RIJ,RIJSQ,A2,DK
 REAL(D):: UIJ

 UIJ=A2*DEXP(-DK*RIJ)/RIJ
 FORCE=UIJ*( 1.0D0/RIJSQ + DK/RIJ)

 RETURN
END FUNCTION FORCE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION BETAURWL(SIGMA12,ESTAR12,DZW)
 IMPLICIT NONE
 INTEGER, PARAMETER:: D  = KIND(1.0D0)
 REAL(D), PARAMETER:: AL = 3.07002D0
 REAL(D):: BETAURWL,SIGMA12,ESTAR12,DZW
 REAL(D):: Z4,Z10

 Z4=(SIGMA12/DZW)**4
 Z10=(SIGMA12/DZW)**10
 BETAURWL=AL*ESTAR12*(Z10 - Z4) + ESTAR12
 
 RETURN
END FUNCTION BETAURWL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION FORCEWL(SIGMA12,ESTAR12,DZW)
 IMPLICIT NONE
 INTEGER, PARAMETER:: D  = KIND(1.0D0)
 REAL(D), PARAMETER:: AL = 3.07002D0
 REAL(D):: FORCEWL,SIGMA12,ESTAR12,DZW
 REAL(D):: IZ,Z4,Z10
 
 IZ=1.0D0/DZW
 Z4=(SIGMA12*IZ)**4
 Z10=(SIGMA12*IZ)**10
 FORCEWL=AL*ESTAR12*(10.0D0*Z10*IZ - 4.0D0*Z4*IZ)

 RETURN
END FUNCTION FORCEWL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


