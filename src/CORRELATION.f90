 INTEGER, PARAMETER:: OUTNSAM = 10           !CORRELATION TIME
 INTEGER, PARAMETER:: TMAXD   = 100000      !MAXIMUM NUMBER OF INSTANT TIMES
 INTEGER, PARAMETER:: IT0     = 10          !TAKE A NEW t=0
 INTEGER, PARAMETER:: T0MAX   = 200000      !MAXIMUM NUMBER OF t=0

 REAL(D), PARAMETER:: TOBS   = 30.0         !OBSERVATION TIME


 REAL(D):: RX0(2,TMAXD),RY0(2,TMAXD),RZ0(2,TMAXD)
 REAL(D):: MSDX(TMAXD,2),MSDY(TMAXD,2),MSDZ(TMAXD,2),MSDT(TMAXD,2)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE CORRELATION
 USE MDVAR
 IMPLICIT NONE
 INTEGER:: I,T,K
 REAL(D):: CTIME,DTIME,TA


 IF(SWITCH .EQ. 0)THEN                      !SET TO ZERO CORRELATION FUNCTIONS
   NTEL=0
   DTIME=H*DBLE(OUTNSAM)                    !CORRELATION TIME STEP
   TMAX=INT(TOBS/DTIME)                     !TOTAL NUMBER OF TIME STEP
   DO I=1,TMAX
      NTIME(I)=0                            !NUMBER OF SAMPLE FOR TIME t_i
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
        Cambiar lado derecho
        RX0(I,TT0)=RX(I)                  !RX(t=0)
        RY0(I,TT0)=RY(K-I)                  !RY(t=0)
        RZ0(I,TT0)=RZ(K-I)                  !RZ(t=0)

     ENDDO    
   ENDIF

   DO T=1,MIN(T0,T0MAX)
      DELT=NTEL - TIME0(T) + 1
      IF(DELT .LE. TMAX)THEN
        NTIME(DELT)=NTIME(DELT) + 1
        DO I=1,2
           K=NPART+1 
           !<[x(0) - x(t)]^2>
           MSDX(DELT,I)=MSDX(DELT,I) + (RX(K-I) - RX0(I,T))**2
           MSDY(DELT,I)=MSDY(DELT,I) + (RY(K-I) - RY0(I,T))**2
           MSDZ(DELT,I)=MSDZ(DELT,I) + (RZ(K-I) - RZ0(I,T))**2
           MSDT(DELT,I)=MSDT(DELT,I) + (RX(K-I) - RX0(I,T))**2 &
                                     + (RY(K-I) - RY0(I,T))**2 &
                                     + (RZ(K-I) - RZ0(I,T))**2
        ENDDO
      ENDIF
   ENDDO
 ELSEIF(SWITCH .EQ. 2)THEN                  !COMPUTE AVERAGES
   OPEN(UNIT=15,FILE='MSD.dat')
   
   DO I=1,TMAX                              !AVERAGE OF CORRELATION FUNCTIONS
      DO K=1,2
         MSDX(I,K)=MSDX(I,K)/DBLE(NTIME(I)*NPART)
         MSDY(I,K)=MSDY(I,K)/DBLE(NTIME(I))
         MSDZ(I,K)=MSDZ(I,K)/DBLE(NTIME(I))
         MSDT(I,K)=MSDT(I,K)/DBLE(NTIME(I))
      ENDDO
      CTIME=H*REAL(OUTNSAM)*(REAL(I))
      WRITE(15,*)CTIME,MSDX(I,K),MSDY(I,K),MSDZ(I,K),MSDT(I,K)
   ENDDO




   CLOSE(UNIT=15)
 ENDIF

 RETURN
END SUBROUTINE CORRELATION
