! Distorsiona una imagen
        SUBROUTINE UNRECTIFY(NEWBUFF)
        IMPLICIT NONE
        INTEGER NEWBUFF
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER NMAPMAX
        PARAMETER (NMAPMAX=3) !second-degree bivariate approximation
        INTEGER NECUAMAX
        PARAMETER (NECUAMAX=(NMAPMAX+1)*NMAPMAX)
        INTEGER MXCOR
        PARAMETER (MXCOR=8) !numero maximo de intersecciones de dos pixels
!
        INTEGER READILIM
!
        INTEGER I,J,II,JJ
        INTEGER I1,I2,J1,J2
        INTEGER NCBUFF
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        INTEGER NMAP
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER IXMIN,IXMAX,IYMIN,IYMAX
        INTEGER NCOR
        INTEGER NEXTINFO
        INTEGER NXMAX_,NYMAX_
        INTEGER IRESAMPL
        REAL AIJ(NECUAMAX),BIJ(NECUAMAX)
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL U,V,X,Y
        REAL U1,V1,U2,V2
        REAL X1,X2,X3,X4
        REAL Y1,Y2,Y3,Y4
        REAL XORIG,YORIG
        REAL XCOR1(2,4),XCOR2(2,4),XCOM(2,MXCOR)
        REAL FAREA,FAREATOT
        CHARACTER*50 CDUMMY
        LOGICAL LERR
!
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN2/NCBUFF
        COMMON/BLKMAPPING0/NMAP
        COMMON/BLKMAPPING1/AIJ,BIJ
        COMMON/BLKMAPPING2/AIJ_,BIJ_
        COMMON/BLKMAPPING3/SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        COMMON/BLKMAPPING4/SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        COMMON/BLKDEFAULTS5/IRESAMPL
!------------------------------------------------------------------------------
        LERR=(NCBUFF.LE.NMAXBUFF/2).AND.(NEWBUFF.LE.NMAXBUFF/2)
!------------------------------------------------------------------------------
        I1=SYMINGRID-SYMINEXTG !borde inferior del primer pixel
        I2=SYMAXGRID+SYMAXEXTG !borde superior del ultimo pixel
        J1=SXMINGRID-SXMINEXTG !borde izquierdo del primer pixel
        J2=SXMAXGRID+SXMAXEXTG !borde derecho del ultimo pixel
!------------------------------------------------------------------------------
        WRITE(CDUMMY,*) NXMAX
        NXMAX_=READILIM('NXMAX of unrectified image',CDUMMY,1,NXMAX)
        WRITE(CDUMMY,*) NYMAX
        NYMAX_=READILIM('NYMAX of unrectified image',CDUMMY,1,NYMAX)
!
        NAXIS(1,NEWBUFF)=NXMAX_
        NAXIS(2,NEWBUFF)=NYMAX_
        DO I=1,NAXIS(2,NEWBUFF)
          DO J=1,NAXIS(1,NEWBUFF)
            IMAGEN(J,I,NEWBUFF)=0.
          END DO
        END DO
        IF(LERR)THEN
          NAXIS(1,NEWBUFF+NMAXBUFF/2)=NXMAX_
          NAXIS(2,NEWBUFF+NMAXBUFF/2)=NYMAX_
          DO I=1,NAXIS(2,NEWBUFF+NMAXBUFF/2)
            DO J=1,NAXIS(1,NEWBUFF+NMAXBUFF/2)
              IMAGEN(J,I,NEWBUFF+NMAXBUFF/2)=0.
            END DO
          END DO
        END IF
!------------------------------------------------------------------------------
        NEXTINFO=0
!------------------------------------------------------------------------------
! Nearest Neighbor
        IF(IRESAMPL.EQ.1)THEN
          DO I=I1,I2-1
            V=REAL(I)+0.5
            DO J=J1,J2-1
              U=REAL(J)+0.5
              CALL FMAP(NMAP,AIJ_,BIJ_,U,V,X,Y)
              JJ=NINT(X)
              II=NINT(Y)
              IMAGEN(JJ,II,NEWBUFF)=IMAGEN(J-J1+1,I-I1+1,NCBUFF)
              IF(LERR)THEN
                IMAGEN(JJ,II,NEWBUFF+NMAXBUFF/2)=IMAGEN(J-J1+1,I-I1+1,NCBUFF+NMAXBUFF/2)
              END IF
            END DO
            CALL SHOWPERC(I1,I2-1,1,I,NEXTINFO)
          END DO
          RETURN
        END IF
!------------------------------------------------------------------------------
! Metodo refinado
        DO I=I1,I2-1
          V1=REAL(I) !borde inferior
          V2=V1+1.0  !borde superior
          DO J=J1,J2-1
            U1=REAL(J) !borde izquierdo
            U2=U1+1.0 !borde derecho
            CALL FMAP(NMAP,AIJ_,BIJ_,U1,V1,X1,Y1)
            CALL FMAP(NMAP,AIJ_,BIJ_,U2,V1,X2,Y2)
            CALL FMAP(NMAP,AIJ_,BIJ_,U2,V2,X3,Y3)
            CALL FMAP(NMAP,AIJ_,BIJ_,U1,V2,X4,Y4)
            XORIG=X1
            YORIG=Y1
            XCOR2(1,1)=0.
            XCOR2(2,1)=0.
            XCOR2(1,2)=X2-XORIG
            XCOR2(2,2)=Y2-YORIG
            XCOR2(1,3)=X3-XORIG
            XCOR2(2,3)=Y3-YORIG
            XCOR2(1,4)=X4-XORIG
            XCOR2(2,4)=Y4-YORIG
            CALL PPOLYN(4,XCOR2,FAREATOT)
            IXMIN=MIN0(NINT(X1),NINT(X2)) !primer pixel a considerar en X
            IXMIN=MIN0(IXMIN,NINT(X3))
            IXMIN=MIN0(IXMIN,NINT(X4))
            IXMAX=MAX0(NINT(X1),NINT(X2)) !ultimo pixel a considerar en X
            IXMAX=MAX0(IXMAX,NINT(X3))
            IXMAX=MAX0(IXMAX,NINT(X4))
            IYMIN=MIN0(NINT(Y1),NINT(Y2)) !primer pixel a considerar en Y
            IYMIN=MIN0(IYMIN,NINT(Y3))
            IYMIN=MIN0(IYMIN,NINT(Y4))
            IYMAX=MAX0(NINT(Y1),NINT(Y2)) !ultimo pixel a considerar en Y
            IYMAX=MAX0(IYMAX,NINT(Y3))
            IYMAX=MAX0(IYMAX,NINT(Y4))
            IF(IXMIN.LT.1) IXMIN=1
            IF(IXMAX.GT.NAXIS(1,NEWBUFF)) IXMAX=NAXIS(1,NEWBUFF)
            IF(IYMIN.LT.1) IYMIN=1
            IF(IYMAX.GT.NAXIS(2,NEWBUFF)) IYMAX=NAXIS(2,NEWBUFF)
            DO II=IYMIN,IYMAX
              DO JJ=IXMIN,IXMAX
                X1=REAL(JJ)-0.5
                X2=X1+1.
                Y1=REAL(II)-0.5
                Y2=Y1+1.
                XCOR1(1,1)=X1-XORIG
                XCOR1(2,1)=Y1-YORIG
                XCOR1(1,2)=X2-XORIG
                XCOR1(2,2)=XCOR1(2,1)
                XCOR1(1,3)=XCOR1(1,2)
                XCOR1(2,3)=Y2-YORIG
                XCOR1(1,4)=XCOR1(1,1)
                XCOR1(2,4)=XCOR1(2,3)
                CALL SRCPAT(4,XCOR1,4,XCOR2,NCOR,XCOM)
                CALL PPOLYN(NCOR,XCOM,FAREA)
                IMAGEN(JJ,II,NEWBUFF)=IMAGEN(JJ,II,NEWBUFF)+IMAGEN(J-J1+1,I-I1+1,NCBUFF)*FAREA/FAREATOT
                IF(LERR)THEN !sumamos los errores al cuadrado
                  IMAGEN(JJ,II,NEWBUFF+NMAXBUFF/2)=IMAGEN(JJ,II,NEWBUFF+NMAXBUFF/2)+ &
                   IMAGEN(J-J1+1,I-I1+1,NCBUFF+NMAXBUFF/2)*IMAGEN(J-J1+1,I-I1+1,NCBUFF+NMAXBUFF/2)*FAREA*FAREA/(FAREATOT*FAREATOT)
                END IF
              END DO
            END DO
          END DO
          CALL SHOWPERC(I1,I2-1,1,I,NEXTINFO)
        END DO
!
        IF(LERR)THEN
          DO I=1,NAXIS(2,NEWBUFF+NMAXBUFF/2)
            DO J=1,NAXIS(1,NEWBUFF+NMAXBUFF/2)
              IMAGEN(J,I,NEWBUFF+NMAXBUFF/2)=SQRT(IMAGEN(J,I,NEWBUFF+NMAXBUFF/2))
            END DO
          END DO
        END IF
!------------------------------------------------------------------------------
        END
