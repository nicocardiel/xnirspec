! Rectifica la imagen actual y la almacena en el buffer NEWBUFF. En caso
! necesario, tambien rectifica la imagen de errores.
        SUBROUTINE RECTIFY(NEWBUFF)
        IMPLICIT NONE
        INTEGER NEWBUFF
!
        INCLUDE 'dimensions.inc'
!
        INTEGER NMAPMAX
        PARAMETER (NMAPMAX=3) !second-degree bivariate approximation
        INTEGER NECUAMAX
        PARAMETER (NECUAMAX=(NMAPMAX+1)*NMAPMAX)
        INTEGER MXCOR
        PARAMETER (MXCOR=8) !numero maximo de intersecciones de dos pixels
!
        INTEGER I,I1,I2,J,J1,J2,II,JJ
        INTEGER IXMIN,IXMAX,IYMIN,IYMAX
        INTEGER NMAP
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER IRESAMPL
        INTEGER NCOR
        INTEGER NEXTINFO
        INTEGER NCBUFF
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL AIJ(NECUAMAX),BIJ(NECUAMAX)
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
        REAL X,Y,U,V
        REAL X1,Y1,X2,Y2,X3,Y3,X4,Y4
        REAL XORIG,YORIG
        REAL U1,V1,U2,V2
        REAL XCOR1(2,4),XCOR2(2,4),XCOM(2,MXCOR)
        REAL FAREA,FAREATOT
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
        IF(NEWBUFF.EQ.NCBUFF)THEN
          WRITE(*,101) '***ERROR***'
          WRITE(*,101) '=> in subroutine RECTIFY'
          WRITE(*,100) '=> NCBUFF,NEWBUFF: '
          WRITE(*,*) NCBUFF,NEWBUFF
          WRITE(*,101) '=> NCBUFF and NEWBUFF must be different!'
          WRITE(*,100) '(press <CR> to continue...)'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
        LERR=(NCBUFF.LE.NMAXBUFF/2).AND.(NEWBUFF.LE.NMAXBUFF/2)
        NEXTINFO=0
!------------------------------------------------------------------------------
        I1=SYMINGRID-SYMINEXTG !borde inferior del primer pixel
        I2=SYMAXGRID+SYMAXEXTG !borde superior del ultimo pixel
        J1=SXMINGRID-SXMINEXTG !borde izquierdo del primer pixel
        J2=SXMAXGRID+SXMAXEXTG !borde derecho del ultimo pixel
! OJO: el numero de pixels no es I2-I1+1, sino que es I2-I1; por eso los bucles
! van desde I1 hasta I2-1. I1,I2,J1,J2 se refieren a los bordes de los pixels,
! no al numero de pixel. Hay un offset de 0.5 entre el mapeado realizado con
! los coeficientes AIJ,BIJ,AIJ_,BIJ_ y el numero de pixel en la imagen original.
        NAXIS(1,NEWBUFF)=J2-J1
        NAXIS(2,NEWBUFF)=I2-I1
        IF(LERR)THEN
          NAXIS(1,NEWBUFF+NMAXBUFF/2)=J2-J1
          NAXIS(2,NEWBUFF+NMAXBUFF/2)=I2-I1
        END IF
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
              IMAGEN(J-J1+1,I-I1+1,NEWBUFF)=IMAGEN(JJ,II,NCBUFF)
              IF(LERR)THEN
                IMAGEN(J-J1+1,I-I1+1,NEWBUFF+NMAXBUFF/2)=IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)
              END IF
            END DO
            CALL SHOWPERC(I1,I2-1,1,I,NEXTINFO)
          END DO
!------------------------------------------------------------------------------
! Linear interpolation: aqui asumimos que los pixels distorsionados son,
! localmente, cuadrilateros (lados rectos). El error va a ser diminuto y de
! hecho tenemos mucha mas incertidumbre con el reparto de la sen~al en el
! pixel.
        ELSEIF(IRESAMPL.EQ.2)THEN
          DO I=I1,I2-1
            V1=REAL(I) !borde inferior
            V2=V1+1.0  !borde superior
            DO J=J1,J2-1
              U1=REAL(J) !borde izquierdo
              U2=U1+1.0  !borde derecho
              CALL FMAP(NMAP,AIJ_,BIJ_,U1,V1,X1,Y1)
              CALL FMAP(NMAP,AIJ_,BIJ_,U2,V1,X2,Y2)
              CALL FMAP(NMAP,AIJ_,BIJ_,U2,V2,X3,Y3)
              CALL FMAP(NMAP,AIJ_,BIJ_,U1,V2,X4,Y4)
!!!              CALL PPOLY4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,3)
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
              IF(IXMAX.GT.NAXIS(1,NCBUFF)) IXMAX=NAXIS(1,NCBUFF)
              IF(IYMIN.LT.1) IYMIN=1
              IF(IYMAX.GT.NAXIS(2,NCBUFF)) IYMAX=NAXIS(2,NCBUFF)
              IMAGEN(J-J1+1,I-I1+1,NEWBUFF)=0.
              IF(LERR) IMAGEN(J-J1+1,I-I1+1,NEWBUFF+NMAXBUFF/2)=0.
              DO II=IYMIN,IYMAX
                DO JJ=IXMIN,IXMAX
                  X1=REAL(JJ)-0.5
                  X2=X1+1.
                  Y1=REAL(II)-0.5
                  Y2=Y1+1.
!!!                  CALL PPOLY4(X1,Y1,X2,Y1,X2,Y2,X1,Y2,6)
                  XCOR1(1,1)=X1-XORIG
                  XCOR1(2,1)=Y1-YORIG
                  XCOR1(1,2)=X2-XORIG
                  XCOR1(2,2)=XCOR1(2,1)
                  XCOR1(1,3)=XCOR1(1,2)
                  XCOR1(2,3)=Y2-YORIG
                  XCOR1(1,4)=XCOR1(1,1)
                  XCOR1(2,4)=XCOR1(2,3)
!!!                  CALL PPOLYN(4,XCOR1,FAREATOT)
                  FAREATOT=1.0 !lo sabemos, por lo que no hacemos PPOLYN
                  CALL SRCPAT(4,XCOR1,4,XCOR2,NCOR,XCOM)
                  IF(NCOR.GT.0)THEN
                    CALL PPOLYN(NCOR,XCOM,FAREA)
                    IMAGEN(J-J1+1,I-I1+1,NEWBUFF)=IMAGEN(J-J1+1,I-I1+1,NEWBUFF)+IMAGEN(JJ,II,NCBUFF)*FAREA/FAREATOT
                    IF(LERR)THEN !sumamos los errores al cuadrado
                      IMAGEN(J-J1+1,I-I1+1,NEWBUFF+NMAXBUFF/2)=IMAGEN(J-J1+1,I-I1+1,NEWBUFF+NMAXBUFF/2)+ &
                       IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)*IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)*FAREA*FAREA/(FAREATOT*FAREATOT)
                    END IF
                  END IF
!!!               read(*,*)
                END DO
              END DO
              IF(LERR)THEN
                IMAGEN(J-J1+1,I-I1+1,NEWBUFF+NMAXBUFF/2)=SQRT(IMAGEN(J-J1+1,I-I1+1,NEWBUFF+NMAXBUFF/2))
              END IF
            END DO
            CALL SHOWPERC(I1,I2-1,1,I,NEXTINFO)
          END DO
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
!
!******************************************************************************
! Dibuja el poligono de 4 lados cuyos vertices se dan como parametros de
! entrada
        SUBROUTINE PPOLY4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,NCOLOR)
        IMPLICIT NONE
        REAL X1,Y1,X2,Y2,X3,Y3,X4,Y4
        INTEGER NCOLOR
!
        INTEGER OLDCI
!------------------------------------------------------------------------------
        CALL PGQCI(OLDCI)
        CALL PGSCI(NCOLOR)
        CALL PGMOVE(X1,Y1)
        CALL PGDRAW(X2,Y2)
        CALL PGDRAW(X3,Y3)
        CALL PGDRAW(X4,Y4)
        CALL PGDRAW(X1,Y1)
        CALL PGSCI(OLDCI)
!
        END
!
!******************************************************************************
! Dibuja la interseccion del pixel problema sobre el pixel original, utilizando
! la informacion calculada por la subrutina SRCPAT. Esta subrutina devuelve
! el area de la interseccion.
        SUBROUTINE PPOLYN(NCOR,XCOM,FAREA)
        IMPLICIT NONE
        INTEGER NCOR
        REAL XCOM(2,NCOR)
        REAL FAREA
!
        INTEGER I1,I2
!!!        INTEGER I,OLDCI
!------------------------------------------------------------------------------
!!!        CALL PGQCI(OLDCI)
!!!        CALL PGSCI(4)
!!!        CALL PGMOVE(XCOM(1,1),XCOM(2,1))
!!!        DO I=2,NCOR
!!!          CALL PGDRAW(XCOM(1,I),XCOM(2,I))
!!!        END DO
!!!        CALL PGDRAW(XCOM(1,1),XCOM(2,1))
!!!        CALL PGSCI(OLDCI)
!
! Nota: es importante tener valores cercanos a uno porque cuando el valor
! del pixel es grande, el valor de FAREA tiene un error muy grande (incluso
! da valores negativos). El modulo que llama a esta subrutina ya se encarga
! de pasar en XCOM valores cercanos a uno.
        FAREA=0.
        DO I1=1,NCOR
          IF(I1.EQ.NCOR)THEN
            I2=1
          ELSE
            I2=I1+1
          END IF
          FAREA=FAREA+XCOM(1,I1)*XCOM(2,I2)-XCOM(1,I2)*XCOM(2,I1)
        END DO
        FAREA=FAREA/2.
!
        END
