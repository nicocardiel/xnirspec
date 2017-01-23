C Dibuja una escala con la calibracion en longitud de onda
        SUBROUTINE PWAVE(NDEGWV,CWV)
        IMPLICIT NONE
        INTEGER NDEGWV
        REAL CWV(20)
C
        INTEGER NMAPMAX
        PARAMETER (NMAPMAX=3) !second-degree bivariate approximation
        INTEGER NECUAMAX
        PARAMETER (NECUAMAX=(NMAPMAX+1)*NMAPMAX)
C
        INCLUDE 'largest.inc'
C
        REAL FPOLY
C
        INTEGER I,L
        INTEGER I1,I2,J1,J2
        INTEGER NSUBX
        INTEGER IDUM,IXINT,IXINT2
        INTEGER NLINEA
        INTEGER NMAP
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        INTEGER NX1,NX2,NY1,NY2
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
        REAL XLEN0,YLEN0
        REAL XMIN,XMAX
        REAL VPX1,VPX2,VPY1,VPY2
        REAL FACTOR
        REAL X0,Y0,DY0
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL XINT,XINT2
        REAL PGRND
        REAL XCH,YCH
        CHARACTER*50 CDUMMY
C
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKMAPPING0/NMAP
        COMMON/BLKMAPPING2/AIJ_,BIJ_
        COMMON/BLKMAPPING3/SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        COMMON/BLKMAPPING4/SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
C------------------------------------------------------------------------------
C limites del mapping
        I1=SYMINGRID-SYMINEXTG !borde inferior del primer pixel
        I2=SYMAXGRID+SYMAXEXTG !borde superior del ultimo pixel
        J1=SXMINGRID-SXMINEXTG !borde izquierdo del primer pixel
        J2=SXMAXGRID+SXMAXEXTG !borde derecho del ultimo pixel
        XMIN=FPOLY(NDEGWV,CWV,REAL(J1)+0.5+REAL(SXMINEXTG)) !limite en l.d.o.
        XMAX=FPOLY(NDEGWV,CWV,REAL(J2)+0.5+REAL(SXMINEXTG)) !limite en l.d.o.
C------------------------------------------------------------------------------
        CALL PGLEN(0,'0',XLEN0,YLEN0) !dimensiones del caracter 0
        CALL PGQVP(0,VPX1,VPX2,VPY1,VPY2) !dimensiones del viewport
        XINT=MAX(0.05,MIN(0.7*XLEN0/(VPX2-VPX1),0.20))*(XMAX-XMIN)
        XINT=PGRND(XINT,NSUBX)
        XINT2=XINT/NSUBX
        FACTOR=10. !afinamos a la decima de Angstrom. Nota: si se cambia este
                   !numero tambien hay que modificar PGNUMB
        DY0=0.07*REAL(I2-I1)
C------------------------------------------------------------------------------
        CALL PGQCS(4,XCH,YCH)
        DO NLINEA=1,3
c eje
          Y0=REAL(I1)+REAL(NLINEA-1)*REAL(I2-I1)/2.
          DO I=1,NXYMAX
            X0=REAL(J1)+REAL(I-1)/REAL(NXYMAX-1)*
     +       REAL(J2-J1)
            CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(I),YP(I))
          END DO
          CALL PGSCI(3)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
c ticks grandes
          IDUM=NINT(XMIN*FACTOR)
          IXINT=NINT(XINT*FACTOR)
          IF(MOD(IDUM,IXINT).NE.0) IDUM=(IDUM/IXINT)*IXINT
          X0=REAL(J1)
          DO WHILE(IDUM.LT.NINT(XMAX*FACTOR))
            IDUM=IDUM+IXINT
            IF(IDUM.LT.NINT(XMAX*FACTOR))THEN
              CALL PGNUMB(IDUM,-1,1,CDUMMY,L)
              CALL FPOLYINV(NDEGWV,CWV,REAL(IDUM)/FACTOR,X0,X0)
              X0=X0-0.5-REAL(SXMINEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0-DY0,XP(1),YP(1))
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0+DY0,XP(2),YP(2))
              CALL PGSCI(3)
              CALL PGMOVE(XP(1),YP(1))
              CALL PGDRAW(XP(2),YP(2))
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0,XP(1),YP(1))
              IF((NINT(XP(1)).GE.NX1).AND.(NINT(XP(1)).LE.NX2).AND.
     +           (NINT(YP(1)).GE.NY1).AND.(NINT(YP(1)).LE.NY2))THEN
                CALL PGSCI(5)
                IF(NLINEA.EQ.1)THEN
                  CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0-1.2*DY0,XP(3),YP(3))
                  CALL PGPTXT(XP(3)+XCH/3.,YP(3),90.,1.0,CDUMMY(1:L))
                ELSEIF(NLINEA.EQ.3)THEN
ccc                  CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0+1.2*DY0,XP(3),YP(3))
ccc                  CALL PGPTXT(XP(3)+XCH/3.,YP(3),90.,0.0,CDUMMY(1:L))
                END IF
              END IF
              CALL PGSCI(1)
            END IF
          END DO
c ticks pequen~os
          IDUM=NINT(XMIN*FACTOR)
          IXINT2=NINT(XINT2*FACTOR)
          IF(MOD(IDUM,IXINT2).NE.0) IDUM=(IDUM/IXINT2)*IXINT2
          X0=REAL(J1)
          DO WHILE(IDUM.LT.NINT(XMAX*FACTOR))
            IDUM=IDUM+IXINT2
            IF(IDUM.LT.NINT(XMAX*FACTOR))THEN
              CALL FPOLYINV(NDEGWV,CWV,REAL(IDUM)/FACTOR,X0,X0)
              X0=X0-0.5-REAL(SXMINEXTG)
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0-0.5*DY0,XP(1),YP(1))
              CALL FMAP(NMAP,AIJ_,BIJ_,X0,Y0+0.5*DY0,XP(2),YP(2))
              CALL PGSCI(3)
              CALL PGMOVE(XP(1),YP(1))
              CALL PGDRAW(XP(2),YP(2))
              CALL PGSCI(1)
            END IF
          END DO
        END DO
C------------------------------------------------------------------------------
        END
