SUBROUTINE VERLET
    USE BDVAR
    IMPLICIT NONE
    INTEGER :: I, J
    REAL(DP) :: DX, DY, DZ, RIJ

    NPAIR = 0

    RV = RCUT + 0.3
    DO I = 1, NPART - 1
        DO J = I + 1, NPART
            DX = RX(I) - RX(J)
            DY = RY(I) - RY(J)
            DZ = RZ(I) - RZ(J)

            DX = DX - BOXX * ANINT(DX * IBOXX)  ! CONDICIONES DE FRONTERA
            DY = DY - BOXY * ANINT(DY * IBOXY)
            DZ = DZ - BOXZ * ANINT(DZ * IBOXZ)

            RIJ = DSQRT(DX*DX + DY*DY + DZ*DZ)
            IF (RIJ .LE. RV) THEN
                NPAIR = NPAIR + 1
                VLI(NPAIR) = I
                VLJ(NPAIR) = J
            END IF
        END DO
    END DO
END SUBROUTINE VERLET
