SUBROUTINE GdR
    IMPLICIT NONE
    USE BDVAR
    INTEGER :: NGR, IG
    REAL(DP) :: DELG

    IF (SWITCH1 .EQ. 0) THEN  ! INICIALIZACION
        NGR = 0
        DELG = BOXX / (2 * NHIS)
        DO I = 1, NHIS
            G(I) = 0.0
        END DO
    END IF

    IF (SWITCH1 .EQ. 1) THEN  ! COMENZAR A TOMAR MUESTRAS
        NGR = NGR + 1
        DO I = 1, NPART - 1
            DO J = I + 1, NPART
                XR = X(I) - X(J)
                YR = Y(I) - Y(J)
                ZR = Z(I) - Z(J)

                XR = XR - ANINT(RXN * IBOXX) * BOXX  ! CONDICIONES DE FRONTERA
                YR = YR - ANINT(RYN * IBOXY) * BOXY
                ZR = ZR - ANINT(RZN * IBOXZ) * BOXZ

                R = DSQRT(XR*XR + YR*YR + ZR*ZR)
                IF (R .LT. (BOXX / 2.0)) THEN  ! RADIO A LA MITAD DE LA CAJA
                    IG = INT(R / DELG)
                    G(IG) = G(IG) + 2.0
                END IF
            END DO
        END DO
    END IF

    IF (SWITCH1 .EQ. 2) THEN  ! CALCULAR GdR (ESTADISTICA)
        DO I = 1, NHIS
            R = DELG * (I + 0.5)
            VBE = ((I + 1.0)**3 - I**3) * DELG**3  ! VOLUMEN DONDE CALCULAREMOS GdR
            NID = (4.0 / 3.0) * PI * VBE * RHOSTAR
            G(I) = G(I) / (NGR * NPART * NID)
        END DO
    END IF
END SUBROUTINE GdR
