C Dibuja la linea en la direccion espectral del mapping que pasa por el punto
C de coordenadas X,Y. La subrutina devuelve de paso las coordenadas U,V.
        SUBROUTINE PXMAP(X,Y,U,V)
        IMPLICIT NONE
        REAL X,Y
        REAL U,V
C
        INCLUDE 'largest.inc'
C
        INTEGER NMAPMAX
        PARAMETER (NMAPMAX=3) !second-degree bivariate approximation
        INTEGER NECUAMAX
        PARAMETER (NECUAMAX=(NMAPMAX+1)*NMAPMAX)
C
        INTEGER I
        INTEGER NMAP
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        REAL AIJ(NECUAMAX),BIJ(NECUAMAX)
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL UU,UMIN,UMAX
C
        COMMON/BLKMAPPING0/NMAP
        COMMON/BLKMAPPING1/AIJ,BIJ
        COMMON/BLKMAPPING2/AIJ_,BIJ_
        COMMON/BLKMAPPING3/SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        COMMON/BLKMAPPING4/SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
C------------------------------------------------------------------------------
        UMIN=REAL(SXMINGRID-SXMINEXTG)
        UMAX=REAL(SXMAXGRID+SXMAXEXTG)
        CALL FMAP(NMAP,AIJ,BIJ,X,Y,U,V)
        DO I=1,NXYMAX
           UU=UMIN+REAL(I-1)/REAL(NXYMAX-1)*(UMAX-UMIN)
           CALL FMAP(NMAP,AIJ_,BIJ_,UU,V,XP(I),YP(I))
        END DO
        CALL PGLINE(NXYMAX,XP,YP)
C------------------------------------------------------------------------------
        END
