! resta el cielo en una imagen con distorsion
        SUBROUTINE SKYSUB(NEWBUFF1,NEWBUFF2)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
! subroutine arguments
        INTEGER NEWBUFF1,NEWBUFF2
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
        REAL FMEAN0,FMEAN1
!!!        REAL FMEAN0E
!
        INTEGER I,J,K,II,JJ
        INTEGER I1,I2,J1,J2
        INTEGER ITEST
        INTEGER IV1,IV2
        INTEGER NCBUFF
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        INTEGER NMAP
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER IXMIN,IXMAX,IYMIN,IYMAX
        INTEGER NCOR
        INTEGER NEXTINFO
        INTEGER NSCANSKY,NSCANOBJ
        INTEGER NOVERSAMP
        INTEGER NITER,NITERTOT
        REAL AIJ(NECUAMAX),BIJ(NECUAMAX)
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL XCUT(NXMAX*NOVERSAMPMAX)
!!!        REAL XCUT_ERR(NXMAX*NOVERSAMPMAX)
        REAL XCUT2(NXMAX*NOVERSAMPMAX),XCUT2_ERR(NXMAX*NOVERSAMPMAX)
        REAL YCUT(NYMAX*NOVERSAMPMAX),YCUT_ERR(NYMAX*NOVERSAMPMAX)
        REAL XC1,YC1,XC2,YC2
        REAL U1,V1,U2,V2
        REAL X1,X2,X3,X4
        REAL Y1,Y2,Y3,Y4
        REAL XORIG,YORIG
        REAL XCOR1(2,4),XCOR2(2,4),XCOM(2,MXCOR)
        REAL XP(NXMAX*NOVERSAMPMAX)
        REAL FAREA,FAREASUBPIXEL
!!!        REAL FAREATOT
        REAL FOVERSAMP
!!!        REAL FSIGMA
        CHARACTER*1 CH
        LOGICAL LERR
        LOGICAL IFSCANSKY(NYMAX)
        LOGICAL IFSCANOBJ(NYMAX)
!
        COMMON/BLKNAXIS/NAXIS
!delete COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN2/NCBUFF
        COMMON/BLKMAPPING0/NMAP
        COMMON/BLKMAPPING1/AIJ,BIJ
        COMMON/BLKMAPPING2/AIJ_,BIJ_
        COMMON/BLKMAPPING3/SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        COMMON/BLKMAPPING4/SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
!------------------------------------------------------------------------------
        IF(NEWBUFF1.EQ.NCBUFF)THEN
          WRITE(*,101) '***ERROR***'
          WRITE(*,101) '=> in subroutine SKYSUB'
          WRITE(*,100) '=> NCBUFF,NEWBUFF1: '
          WRITE(*,*) NCBUFF,NEWBUFF1
          WRITE(*,101) '=> NCBUFF and NEWBUFF1 must be different!'
          WRITE(*,100) '(press <CR> to continue...)'
          READ(*,*)
          RETURN
        END IF
        IF(NEWBUFF2.EQ.NCBUFF)THEN
          WRITE(*,101) '***ERROR***'
          WRITE(*,101) '=> in subroutine SKYSUB'
          WRITE(*,100) '=> NCBUFF,NEWBUFF2: '
          WRITE(*,*) NCBUFF,NEWBUFF2
          WRITE(*,101) '=> NCBUFF and NEWBUFF2 must be different!'
          WRITE(*,100) '(press <CR> to continue...)'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
        LERR=(NCBUFF.LE.NMAXBUFF/2).AND.(NEWBUFF1.LE.NMAXBUFF/2).AND.(NEWBUFF2.LE.NMAXBUFF/2)
!------------------------------------------------------------------------------
        I1=SYMINGRID-SYMINEXTG !borde inferior del primer pixel
        I2=SYMAXGRID+SYMAXEXTG !borde superior del ultimo pixel
        J1=SXMINGRID-SXMINEXTG !borde izquierdo del primer pixel
        J2=SXMAXGRID+SXMAXEXTG !borde derecho del ultimo pixel
! OJO: el numero de pixels no es I2-I1+1, sino que es I2-I1; por eso los bucles
! van desde I1 hasta I2-1. I1,I2,J1,J2 se refieren a los bordes de los pixels,
! no al numero de pixel. Hay un offset de 0.5 entre el mapeado realizado con
! los coeficientes AIJ,BIJ,AIJ_,BIJ_ y el numero de pixel en la imagen original.
!------------------------------------------------------------------------------
! IFSCANSKY indica que espectros de la imagen corregida vamos a utilizar para
! calcular el cielo
        DO I=I1,I2-1
          IFSCANSKY(I-I1+1)=.FALSE.
        END DO
!
        WRITE(*,101) 'Select sky region(s)...'
        CH='A'
        CALL PGSCI(3)
        DO WHILE(CH.NE.'X')
          CALL RPGBAND(0,0,0.,0.,XC1,YC1,CH)
          IF(CH.NE.'X')THEN
            CALL PXMAP(XC1,YC1,U1,V1)
            WRITE(*,*) XC1, YC1
            CALL RPGBAND(0,0,0.,0.,XC2,YC2,CH)
            IF(CH.NE.'X')THEN
              CALL PXMAP(XC2,YC2,U2,V2)
              WRITE(*,*) XC2, YC2
              IF(V2.GT.V1)THEN
                IV1=NINT(V1)
                IV2=NINT(V2)
              ELSE
                IV1=NINT(V2)
                IV2=NINT(V1)
              END IF
              DO I=IV1,IV2
                IF((I.GE.I1).AND.(I.LE.I2-1))THEN
                  IFSCANSKY(I-I1+1)=.TRUE.
                END IF
              END DO
              CALL FMAP(NMAP,AIJ_,BIJ_,REAL(J1),REAL(IV1),XC1,YC1)
              CALL FMAP(NMAP,AIJ_,BIJ_,REAL(J1),REAL(IV2),XC2,YC2)
              CALL PGMOVE(XC1,YC1)
              CALL PGDRAW(XC2,YC2)
              CALL FMAP(NMAP,AIJ_,BIJ_,REAL(J2),REAL(IV1),XC1,YC1)
              CALL FMAP(NMAP,AIJ_,BIJ_,REAL(J2),REAL(IV2),XC2,YC2)
              CALL PGMOVE(XC1,YC1)
              CALL PGDRAW(XC2,YC2)
            END IF
          END IF
        END DO
        WRITE(*,101) '...OK!'
        CALL PGSCI(1)
        K=0
        DO I=I1,I2-1
          II=I-I1+1
          IF(IFSCANSKY(II))THEN
            K=K+1
          END IF
        END DO
        IF(K.EQ.0) RETURN
        NSCANSKY=K
        WRITE(*,100) '=> Number of sky spectra: '
        WRITE(*,*) NSCANSKY
!------------------------------------------------------------------------------
! IFSCANOBJ indica que espectros de la imagen contienen el/los objeto/s de
! interes
        DO I=I1,I2-1
          IFSCANOBJ(I-I1+1)=.FALSE.
        END DO
!
        WRITE(*,101) 'Select object region(s)...'
        CH='A'
        CALL PGSCI(4)
        DO WHILE(CH.NE.'X')
          CALL RPGBAND(0,0,0.,0.,XC1,YC1,CH)
          IF(CH.NE.'X')THEN
            CALL PXMAP(XC1,YC1,U1,V1)
            WRITE(*,*) XC1, YC1
            CALL RPGBAND(0,0,0.,0.,XC2,YC2,CH)
            IF(CH.NE.'X')THEN
              CALL PXMAP(XC2,YC2,U2,V2)
              WRITE(*,*) XC2, YC2
              IF(V2.GT.V1)THEN
                IV1=NINT(V1)
                IV2=NINT(V2)
              ELSE
                IV1=NINT(V2)
                IV2=NINT(V1)
              END IF
              DO I=IV1,IV2
                IF((I.GE.I1).AND.(I.LE.I2-1))THEN
                  IFSCANOBJ(I-I1+1)=.TRUE.
                END IF
              END DO
              CALL FMAP(NMAP,AIJ_,BIJ_,REAL(J1),REAL(IV1),XC1,YC1)
              CALL FMAP(NMAP,AIJ_,BIJ_,REAL(J1),REAL(IV2),XC2,YC2)
              CALL PGMOVE(XC1,YC1)
              CALL PGDRAW(XC2,YC2)
              CALL FMAP(NMAP,AIJ_,BIJ_,REAL(J2),REAL(IV1),XC1,YC1)
              CALL FMAP(NMAP,AIJ_,BIJ_,REAL(J2),REAL(IV2),XC2,YC2)
              CALL PGMOVE(XC1,YC1)
              CALL PGDRAW(XC2,YC2)
            END IF
          END IF
        END DO
        WRITE(*,101) '...OK!'
        CALL PGSCI(1)
        K=0
        DO I=I1,I2-1
          II=I-I1+1
          IF(IFSCANOBJ(II))THEN
            K=K+1
          END IF
        END DO
        IF(K.EQ.0) RETURN
        NSCANOBJ=K
        WRITE(*,100) '=> Number of object spectra: '
        WRITE(*,*) NSCANOBJ
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Extraemos el espectro del cielo, usando sobremuestreado si se solicita.
        NOVERSAMP=READILIM('Oversampling per pixel','1',1,NOVERSAMPMAX)
        FOVERSAMP=1./REAL(NOVERSAMP) !subpixel size (pixel units)
        NITERTOT=READILIM('Number of iterations','1',1,50)
        NITER=1
!------------------------------------------------------------------------------
        WRITE(*,101) 'Extracting initial sky spectrum...'
        NEXTINFO=0
        DO J=1,(J2-J1)*NOVERSAMP
          U1=REAL(J1)+REAL(J-1)*FOVERSAMP !borde izquierdo
          U2=U1+FOVERSAMP                 !borde derecho
          K=0
          DO I=1,(I2-I1)*NOVERSAMP
            ITEST=INT((I-1)/NOVERSAMP)+1
            IF(IFSCANSKY(ITEST))THEN
              V1=REAL(I1)+REAL(I-1)*FOVERSAMP !borde inferior
              V2=V1+FOVERSAMP                 !borde superior
              K=K+1
              YCUT(K)=0.0
!              YCUT_ERR(K)=0.0
              CALL FMAP(NMAP,AIJ_,BIJ_,U1,V1,X1,Y1)
              CALL FMAP(NMAP,AIJ_,BIJ_,U2,V1,X2,Y2)
              CALL FMAP(NMAP,AIJ_,BIJ_,U2,V2,X3,Y3)
              CALL FMAP(NMAP,AIJ_,BIJ_,U1,V2,X4,Y4)
!!!           CALL PPOLY4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,3)
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
              DO II=IYMIN,IYMAX
                DO JJ=IXMIN,IXMAX
                  X1=REAL(JJ)-0.5
                  X2=X1+1.
                  Y1=REAL(II)-0.5
                  Y2=Y1+1.
!!!               CALL PPOLY4(X1,Y1,X2,Y1,X2,Y2,X1,Y2,6)
                  XCOR1(1,1)=X1-XORIG
                  XCOR1(2,1)=Y1-YORIG
                  XCOR1(1,2)=X2-XORIG
                  XCOR1(2,2)=XCOR1(2,1)
                  XCOR1(1,3)=XCOR1(1,2)
                  XCOR1(2,3)=Y2-YORIG
                  XCOR1(1,4)=XCOR1(1,1)
                  XCOR1(2,4)=XCOR1(2,3)
!!!               CALL PPOLYN(4,XCOR1,FAREATOT)
!!!               FAREATOT=1.0 !lo sabemos, por lo que no hacemos PPOLYN
                  CALL SRCPAT(4,XCOR1,4,XCOR2,NCOR,XCOM)
                  IF(NCOR.GT.0)THEN
                    CALL PPOLYN(NCOR,XCOM,FAREA)
                    YCUT(K)=YCUT(K)+IMAGEN(JJ,II,NCBUFF)*FAREA!!!/FAREATOT
!                    IF(LERR)THEN !sumamos los errores al cuadrado
!                      YCUT_ERR(K)=YCUT_ERR(K)+
!     +                 IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)*
!     +                 IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)*
!     +                 FAREA*FAREA!!!/(FAREATOT*FAREATOT)
!                    END IF
                  END IF
                END DO
              END DO
!              IF(LERR)THEN
!                YCUT_ERR(K)=SQRT(YCUT_ERR(K))
!              END IF
            END IF
          END DO ! DO I=1,(I2-I1)*NOVERSAMP
!          IF(LERR)THEN
!            XCUT(J)=FMEAN0E(K,YCUT,YCUT_ERR,FSIGMA,XCUT_ERR(J))
!          ELSE
            XCUT(J)=FMEAN1(K,YCUT)
!            XCUT_ERR(J)=0.0
!          END IF
          CALL SHOWPERC(1,(J2-J1)*NOVERSAMP,1,J,NEXTINFO)
        END DO ! DO J=1,(J2-J1)*NOVERSAMP
!------------------------------------------------------------------------------
! La coordenada X de los nuevos subpixels se refieren al centro de dichos
! subpixels. Por esa razón añadimos el factor 0.5*FOVERSAMP.
        DO J=1,(J2-J1)*NOVERSAMP
          XP(J)=REAL(J1)+REAL(J-1)*FOVERSAMP+0.5*FOVERSAMP
          write(77,*) xp(j),xcut(j)
        END DO
        CALL SUBPLOT((J2-J1)*NOVERSAMP,1,(J2-J1)*NOVERSAMP,XP,XCUT,XP,XCUT,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
         'wavelength direction','signal','averaged sky spectrum',3,201,1.0)
!------------------------------------------------------------------------------
! Creamos una imagen distorsionada usando el espectro de cielo promedio
! anterior.
        NAXIS(1,NEWBUFF1)=NAXIS(1,NCBUFF)
        NAXIS(2,NEWBUFF1)=NAXIS(2,NCBUFF)
        NAXIS(1,NEWBUFF2)=NAXIS(1,NCBUFF)
        NAXIS(2,NEWBUFF2)=NAXIS(2,NCBUFF)
10      CONTINUE
        DO I=1,NAXIS(2,NEWBUFF1)
          DO J=1,NAXIS(1,NEWBUFF1)
            IMAGEN(J,I,NEWBUFF1)=0.
            IMAGEN(J,I,NEWBUFF2)=0.
          END DO
        END DO
!        IF(LERR)THEN
!          NAXIS(1,NEWBUFF1+NMAXBUFF/2)=NAXIS(1,NCBUFF+NMAXBUFF/2)
!          NAXIS(2,NEWBUFF1+NMAXBUFF/2)=NAXIS(2,NCBUFF+NMAXBUFF/2)
!          DO I=1,NAXIS(2,NEWBUFF1+NMAXBUFF/2)
!            DO J=1,NAXIS(1,NEWBUFF1+NMAXBUFF/2)
!              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=0.
!            END DO
!          END DO
!        END IF
!
        WRITE(*,100) 'Iteration #'
        WRITE(*,*) NITER
        WRITE(*,101) 'Creating sky image in sky region...'
        NEXTINFO=0
        DO J=1,(J2-J1)*NOVERSAMP
          U1=REAL(J1)+REAL(J-1)*FOVERSAMP !borde izquierdo
          U2=U1+FOVERSAMP                 !borde derecho
          DO I=1,(I2-I1)*NOVERSAMP
            ITEST=INT((I-1)/NOVERSAMP)+1
            IF(IFSCANSKY(ITEST))THEN
              V1=REAL(I1)+REAL(I-1)*FOVERSAMP !borde inferior
              V2=V1+FOVERSAMP                 !borde superior
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
              CALL PPOLYN(4,XCOR2,FAREASUBPIXEL)
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
              IF(IXMAX.GT.NAXIS(1,NEWBUFF1)) IXMAX=NAXIS(1,NEWBUFF1)
              IF(IYMIN.LT.1) IYMIN=1
              IF(IYMAX.GT.NAXIS(2,NEWBUFF1)) IYMAX=NAXIS(2,NEWBUFF1)
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
                  IMAGEN(JJ,II,NEWBUFF1)=IMAGEN(JJ,II,NEWBUFF1)+XCUT(J)*FAREA/FAREASUBPIXEL
                  IMAGEN(JJ,II,NEWBUFF2)=IMAGEN(JJ,II,NEWBUFF2)+FAREA
!                  IF(LERR)THEN !sumamos los errores al cuadrado
!                    IMAGEN(JJ,II,NEWBUFF1+NMAXBUFF/2)=
!     +               IMAGEN(JJ,II,NEWBUFF1+NMAXBUFF/2)+
!     +               XCUT_ERR(J)*XCUT_ERR(J)*FAREA/FAREASUBPIXEL
!                  END IF
                END DO
              END DO
            END IF
          END DO
          CALL SHOWPERC(1,(J2-J1)*NOVERSAMP,1,J,NEXTINFO)
        END DO
!------------------------------------------------------------------------------
! Extraemos de nuevo un espectro de cielo, pero esta vez usamos el espectro
! inicial para repartir la señal dentro de cada píxel.
        WRITE(*,100) 'Iteration #'
        WRITE(*,*) NITER
        WRITE(*,101) 'Extracting refined sky spectrum...'
        NEXTINFO=0
        DO J=1,(J2-J1)*NOVERSAMP
          U1=REAL(J1)+REAL(J-1)*FOVERSAMP !borde izquierdo
          U2=U1+FOVERSAMP                 !borde derecho
          K=0
          DO I=1,(I2-I1)*NOVERSAMP
            ITEST=INT((I-1)/NOVERSAMP)+1
            IF(IFSCANSKY(ITEST))THEN
              V1=REAL(I1)+REAL(I-1)*FOVERSAMP !borde inferior
              V2=V1+FOVERSAMP                 !borde superior
              K=K+1
              YCUT(K)=0.0
              YCUT_ERR(K)=0.0
              CALL FMAP(NMAP,AIJ_,BIJ_,U1,V1,X1,Y1)
              CALL FMAP(NMAP,AIJ_,BIJ_,U2,V1,X2,Y2)
              CALL FMAP(NMAP,AIJ_,BIJ_,U2,V2,X3,Y3)
              CALL FMAP(NMAP,AIJ_,BIJ_,U1,V2,X4,Y4)
!!!           CALL PPOLY4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,3)
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
              CALL PPOLYN(4,XCOR2,FAREASUBPIXEL)
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
              DO II=IYMIN,IYMAX
                DO JJ=IXMIN,IXMAX
                  X1=REAL(JJ)-0.5
                  X2=X1+1.
                  Y1=REAL(II)-0.5
                  Y2=Y1+1.
!!!               CALL PPOLY4(X1,Y1,X2,Y1,X2,Y2,X1,Y2,6)
                  XCOR1(1,1)=X1-XORIG
                  XCOR1(2,1)=Y1-YORIG
                  XCOR1(1,2)=X2-XORIG
                  XCOR1(2,2)=XCOR1(2,1)
                  XCOR1(1,3)=XCOR1(1,2)
                  XCOR1(2,3)=Y2-YORIG
                  XCOR1(1,4)=XCOR1(1,1)
                  XCOR1(2,4)=XCOR1(2,3)
!!!               CALL PPOLYN(4,XCOR1,FAREATOT)
!!!               FAREATOT=1.0 !lo sabemos, por lo que no hacemos PPOLYN
                  CALL SRCPAT(4,XCOR1,4,XCOR2,NCOR,XCOM)
                  IF(NCOR.GT.0)THEN
                    CALL PPOLYN(NCOR,XCOM,FAREA)
                    !!!la siguiente línea está comentada porque hacía
                    !que el método sólo funcionara cuando la señal del
                    !cielo es positiva; sin embargo, nos interesa poder
                    !utilizarlo cuando quiero restar residuos de cielo,
                    !que pueden tener valores positivos y negativos.
                    !IF(IMAGEN(JJ,II,NEWBUFF1).GT.0)THEN
                      !La clave está en redistribuir la señal en el pixel
                      !IMAGEN(JJ,II,NCBUFF) usando la informacion del espectro
                      !anterior XCUT(J). Para ello calculamos qué fracción
                      !de señal del espectro cabe esperar en el subpixel
                      !considerado, teniendo en cuenta que en la reconstrucción
                      !de la imagen de cielo dicho subpixel contribuyó con una
                      !señal FAREA/FAREASUBPIXEL*XCUT(J) a la señal del pixel
                      !(JJ,II).
                      !También es esencial corregir de efecto de borde, usando
                      !para ello la información almacenada en NEWBUFF2.
                      !
                      !Opción a) no funciona del todo mal... pero hay que
                      !iterar bastante
!                      YCUT(K)=YCUT(K)+XCUT(J)*FAREA/FAREASUBPIXEL+
!     +                 (IMAGEN(JJ,II,NCBUFF)-IMAGEN(JJ,II,NEWBUFF1))*
!     +                 FAREA
                      !
                      !Opción b) intentando mejorar... pero no funciona
!                      YCUT(K)=YCUT(K)+XCUT(J)*FAREA/FAREASUBPIXEL+
!     +                 (IMAGEN(JJ,II,NCBUFF)-IMAGEN(JJ,II,NEWBUFF1))*
!     +                 FAREA/FAREASUBPIXEL*XCUT(J)/IMAGEN(JJ,II,NCBUFF)
                      !
                      !Opción c) otra prueba: algo hemos ganado
                      YCUT(K)=YCUT(K)+XCUT(J)*FAREA/FAREASUBPIXEL+(IMAGEN(JJ,II,NCBUFF)*IMAGEN(JJ,II,NEWBUFF2)- &
                       IMAGEN(JJ,II,NEWBUFF1))*FAREA
                      !Opción z) lo original
!                      YCUT(K)=YCUT(K)+
!     +                 IMAGEN(JJ,II,NCBUFF)*FAREA!!!/FAREATOT
!!!                      IF(LERR)THEN !sumamos los errores al cuadrado
!!!                        YCUT_ERR(K)=YCUT_ERR(K)+
!!!     +                   IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)*
!!!     +                   IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)*
!!!     +                   FAREA*FAREA!!!/(FAREATOT*FAREATOT)
!!!                      END IF
                    !!!la siguiente línea está comentada por la
                    !explicación que se da más arriba
                    !END IF
                  END IF
                END DO
              END DO
!!!              IF(LERR)THEN
!!!                YCUT_ERR(K)=SQRT(YCUT_ERR(K))
!!!              END IF
            END IF
          END DO ! DO I=1,(I2-I1)*NOVERSAMP
          IF(LERR)THEN
            XCUT2(J)=FMEAN0(K,YCUT,XCUT2_ERR(J))
          ELSE
            XCUT2(J)=FMEAN1(K,YCUT)
            XCUT2_ERR(J)=0.0
          END IF
          CALL SHOWPERC(1,(J2-J1)*NOVERSAMP,1,J,NEXTINFO)
        END DO ! DO J=1,(J2-J1)*NOVERSAMP
!------------------------------------------------------------------------------
! La coordenada X de los nuevos subpixels se refieren al centro de dichos
! subpixels. Por esa razón añadimos el factor 0.5*FOVERSAMP.
        DO J=1,(J2-J1)*NOVERSAMP
!          XP(J)=REAL(J1)+REAL(J-1)*FOVERSAMP+0.5*FOVERSAMP
          write(77+NITER,*) xp(j),xcut2(j)
        END DO
        CALL SUBPLOT((J2-J1)*NOVERSAMP,1,(J2-J1)*NOVERSAMP,XP,XCUT2,XP,XCUT2,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
         'wavelength direction','signal','averaged sky spectrum',3,201,1.0)
!------------------------------------------------------------------------------
        IF(NITER.LT.NITERTOT)THEN
          NITER=NITER+1
          DO J=1,(J2-J1)*NOVERSAMP
            XCUT(J)=XCUT2(J)
          END DO
          GOTO 10
        END IF
!------------------------------------------------------------------------------
! Creamos una imagen distorsionada usando el espectro de cielo promedio
! anterior.
        NAXIS(1,NEWBUFF1)=NAXIS(1,NCBUFF)
        NAXIS(2,NEWBUFF1)=NAXIS(2,NCBUFF)
        DO I=1,NAXIS(2,NEWBUFF1)
          DO J=1,NAXIS(1,NEWBUFF1)
            IMAGEN(J,I,NEWBUFF1)=0.
          END DO
        END DO
        IF(LERR)THEN
          NAXIS(1,NEWBUFF1+NMAXBUFF/2)=NAXIS(1,NCBUFF+NMAXBUFF/2)
          NAXIS(2,NEWBUFF1+NMAXBUFF/2)=NAXIS(2,NCBUFF+NMAXBUFF/2)
          DO I=1,NAXIS(2,NEWBUFF1+NMAXBUFF/2)
            DO J=1,NAXIS(1,NEWBUFF1+NMAXBUFF/2)
              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=0.
            END DO
          END DO
        END IF
!
        WRITE(*,101) 'Creating sky image in object region...'
        NEXTINFO=0
        DO J=1,(J2-J1)*NOVERSAMP
          U1=REAL(J1)+REAL(J-1)*FOVERSAMP !borde izquierdo
          U2=U1+FOVERSAMP                 !borde derecho
          DO I=1,(I2-I1)*NOVERSAMP
            ITEST=INT((I-1)/NOVERSAMP)+1
            IF(IFSCANOBJ(ITEST))THEN
              V1=REAL(I1)+REAL(I-1)*FOVERSAMP !borde inferior
              V2=V1+FOVERSAMP                 !borde superior
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
              CALL PPOLYN(4,XCOR2,FAREASUBPIXEL)
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
              IF(IXMAX.GT.NAXIS(1,NEWBUFF1)) IXMAX=NAXIS(1,NEWBUFF1)
              IF(IYMIN.LT.1) IYMIN=1
              IF(IYMAX.GT.NAXIS(2,NEWBUFF1)) IYMAX=NAXIS(2,NEWBUFF1)
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
                  IMAGEN(JJ,II,NEWBUFF1)=IMAGEN(JJ,II,NEWBUFF1)+XCUT2(J)*FAREA/FAREASUBPIXEL
                  IF(LERR)THEN !sumamos los errores al cuadrado
                    IMAGEN(JJ,II,NEWBUFF1+NMAXBUFF/2)=IMAGEN(JJ,II,NEWBUFF1+NMAXBUFF/2)+ &
                     XCUT2_ERR(J)*XCUT2_ERR(J)*FAREA/FAREASUBPIXEL
                  END IF
                END DO
              END DO
            END IF
          END DO
          CALL SHOWPERC(1,(J2-J1)*NOVERSAMP,1,J,NEXTINFO)
        END DO
!
        IF(LERR)THEN
          DO I=1,NAXIS(2,NEWBUFF1+NMAXBUFF/2)
            DO J=1,NAXIS(1,NEWBUFF1+NMAXBUFF/2)
              IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)=SQRT(IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2))
            END DO
          END DO
        END IF
!------------------------------------------------------------------------------
! restamos por fin la imagen
        DO I=1,NAXIS(2,NCBUFF)
          DO J=1,NAXIS(1,NCBUFF)
            IMAGEN(J,I,NCBUFF)=IMAGEN(J,I,NCBUFF)-IMAGEN(J,I,NEWBUFF1)
          END DO
        END DO
        IF(LERR)THEN
          DO I=1,NAXIS(2,NCBUFF+NMAXBUFF/2)
            DO J=1,NAXIS(1,NCBUFF+NMAXBUFF/2)
              IMAGEN(J,I,NCBUFF+NMAXBUFF/2)=SQRT( &
               IMAGEN(J,I,NCBUFF+NMAXBUFF/2)*IMAGEN(J,I,NCBUFF+NMAXBUFF/2)+ &
               IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2)*IMAGEN(J,I,NEWBUFF1+NMAXBUFF/2))
            END DO
          END DO
        END IF
        CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
