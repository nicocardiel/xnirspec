! Retorna el valor en la columna numero NCOL de la tira de numeros en CLINEA
        REAL FUNCTION LEECOLUMN(CLINEA,NCOL)
        IMPLICIT NONE
        CHARACTER*(*) CLINEA
        INTEGER NCOL
!
        INTEGER TRUEBEG,TRUELEN
!
        INTEGER N
        INTEGER L1,L2,L
!------------------------------------------------------------------------------
        L1=TRUEBEG(CLINEA)
        L2=TRUELEN(CLINEA)
        IF((L1.EQ.0).OR.(L2.EQ.0))THEN
          WRITE(*,101) 'ERROR in LEECOLUMN: empty string!'
          LEECOLUMN=0.0
          RETURN
        END IF
! Caso trivial
        IF(NCOL.EQ.1)THEN
          READ(CLINEA,*) LEECOLUMN
          RETURN
        END IF
!
        N=1 !numero de columnas ignoradas
        DO WHILE(N.LT.NCOL)
          L=INDEX(CLINEA(L1:L2),' ') !siguiente espacio en blanco
          IF(L.EQ.0)THEN
            WRITE(*,101)'ERROR in LEECOLUMN: column no. does not exist.'
            RETURN
          END IF
          L=L1+L-1
          L1=TRUEBEG(CLINEA(L:L2))
          L1=L1+L-1
          N=N+1
          IF(N.EQ.NCOL)THEN
            READ(CLINEA(L1:L2),*) LEECOLUMN
            RETURN
          END IF
        END DO
!
101     FORMAT(A)
        END
