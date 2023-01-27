        SUBROUTINE PMAPP
        IMPLICIT NONE
!
        INCLUDE 'largest.inc'
!
!
        INTEGER NMAPMAX
        PARAMETER (NMAPMAX=3) !second-degree bivariate approximation
        INTEGER NECUAMAX
        PARAMETER (NECUAMAX=(NMAPMAX+1)*NMAPMAX)
!
        INTEGER I
        INTEGER NX0,NX1,NY0,NY1
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        INTEGER NGRID,NMAP
        REAL X0,Y0
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL AIJ(NECUAMAX),BIJ(NECUAMAX)
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
!
        COMMON/BLKDEFAULTS4/NGRID
        COMMON/BLKMAPPING0/NMAP
        COMMON/BLKMAPPING1/AIJ,BIJ
        COMMON/BLKMAPPING2/AIJ_,BIJ_
        COMMON/BLKMAPPING3/SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        COMMON/BLKMAPPING4/SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
!------------------------------------------------------------------------------
        CALL PGBBUF
!------------------------------------------------------------------------------
! dibujamos extrapolacion hacia abajo
        IF(SYMINEXTG.GT.0)THEN
          NY0=SYMINGRID
          NY1=SYMINGRID-SYMINEXTG
          DO WHILE(NY0.GE.NY1)
            Y0=REAL(NY0)
            DO I=1,NXYMAX
              X0=REAL(SXMINGRID-SXMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SXMAXGRID+SXMAXEXTG-SXMINGRID+SXMINEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
            END DO
            CALL PGSCI(5)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            NY0=NY0-NGRID
          END DO
          NX0=SXMINGRID
          NX1=SXMAXGRID
          DO WHILE(NX0.LE.NX1)
            X0=REAL(NX0)
            DO I=1,NXYMAX
              Y0=REAL(SYMINGRID-SYMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SYMINEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
            END DO
            CALL PGSCI(5)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            NX0=NX0+NGRID
          END DO
        END IF
! dibujamos extrapolacion hacia arriba
        IF(SYMAXEXTG.GT.0)THEN
          NY0=SYMAXGRID
          NY1=SYMAXGRID+SYMAXEXTG
          DO WHILE(NY0.LE.NY1)
            Y0=REAL(NY0)
            DO I=1,NXYMAX
              X0=REAL(SXMINGRID-SXMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SXMAXGRID+SXMAXEXTG-SXMINGRID+SXMINEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
            END DO
            CALL PGSCI(5)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            NY0=NY0+NGRID
          END DO
          NX0=SXMINGRID
          NX1=SXMAXGRID
          DO WHILE(NX0.LE.NX1)
            X0=REAL(NX0)
            DO I=1,NXYMAX
              Y0=REAL(SYMAXGRID)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SYMAXEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
            END DO
            CALL PGSCI(5)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            NX0=NX0+NGRID
          END DO
        END IF
! dibujamos extrapolacion hacia la izquierda
        IF(SXMINEXTG.GT.0)THEN
          NX0=SXMINGRID
          NX1=SXMINGRID-SXMINEXTG
          DO WHILE(NX0.GE.NX1)
            X0=REAL(NX0)
            DO I=1,NXYMAX
              Y0=REAL(SYMINGRID-SYMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SYMAXGRID+SYMAXEXTG-SYMINGRID+SYMINEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
            END DO
            CALL PGSCI(5)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            NX0=NX0-NGRID
          END DO
          NY0=SYMINGRID
          NY1=SYMAXGRID
          DO WHILE(NY0.LE.NY1)
            Y0=REAL(NY0)
            DO I=1,NXYMAX
              X0=REAL(SXMINGRID-SXMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SXMINEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
            END DO
            CALL PGSCI(5)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            NY0=NY0+NGRID
          END DO
        END IF
! dibujamos extrapolacion hacia la derecha
        IF(SXMAXEXTG.GT.0)THEN
          NX0=SXMAXGRID
          NX1=SXMAXGRID+SXMAXEXTG
          DO WHILE(NX0.LE.NX1)
            X0=REAL(NX0)
            DO I=1,NXYMAX
              Y0=REAL(SYMINGRID-SYMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SYMAXGRID+SYMAXEXTG-SYMINGRID+SYMINEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
            END DO
            CALL PGSCI(5)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            NX0=NX0+NGRID
          END DO
          NY0=SYMINGRID
          NY1=SYMAXGRID
          DO WHILE(NY0.LE.NY1)
            Y0=REAL(NY0)
            DO I=1,NXYMAX
              X0=REAL(SXMAXGRID)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SXMAXEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
            END DO
            CALL PGSCI(5)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            NY0=NY0+NGRID
          END DO
        END IF
!------------------------------------------------------------------------------
! dibujamos una malla equidistante de paso NGRID (usamos variables enteras
! para evitar errores de redondeo al avanzar en el grid)
        NY0=SYMINGRID
        NY1=SYMAXGRID
        DO WHILE(NY0.LE.NY1)
          Y0=REAL(NY0)
          DO I=1,NXYMAX
            X0=REAL(SXMINGRID)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SXMAXGRID-SXMINGRID)
            CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
          END DO
          CALL PGSCI(3)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
          NY0=NY0+NGRID
        END DO
        NX0=SXMINGRID
        NX1=SXMAXGRID
        DO WHILE(NX0.LE.NX1)
          X0=REAL(NX0)
          DO I=1,NXYMAX
            Y0=REAL(SYMINGRID)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SYMAXGRID-SYMINGRID)
            CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
          END DO
          CALL PGSCI(3)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
          NX0=NX0+NGRID
        END DO
!------------------------------------------------------------------------------
! dibujamos perimetro del grid (incluyendo las extrapolaciones)
! borde inferior
          Y0=REAL(SYMINGRID-SYMINEXTG)
          DO I=1,NXYMAX
            X0=REAL(SXMINGRID-SXMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SXMAXGRID+SXMAXEXTG-SXMINGRID+SXMINEXTG)
            CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
          END DO
          CALL PGSCI(7)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
! borde superior
          Y0=REAL(SYMAXGRID+SYMAXEXTG)
          DO I=1,NXYMAX
            X0=REAL(SXMINGRID-SXMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SXMAXGRID+SXMAXEXTG-SXMINGRID+SXMINEXTG)
            CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
          END DO
          CALL PGSCI(7)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
! borde izquierdo
          X0=REAL(SXMINGRID-SXMINEXTG)
          DO I=1,NXYMAX
            Y0=REAL(SYMINGRID-SYMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SYMAXGRID+SYMAXEXTG-SYMINGRID+SYMINEXTG)
            CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
          END DO
          CALL PGSCI(7)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
! borde derecho
          X0=REAL(SXMAXGRID+SXMAXEXTG)
          DO I=1,NXYMAX
            Y0=REAL(SYMINGRID-SYMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SYMAXGRID+SYMAXEXTG-SYMINGRID+SYMINEXTG)
            CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
          END DO
          CALL PGSCI(7)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
!------------------------------------------------------------------------------
        CALL PGEBUF
        END
