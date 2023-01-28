! Extrae un espectro
        SUBROUTINE XEXTRACT
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER NMAPMAX
        PARAMETER (NMAPMAX=3)
        INTEGER NECUAMAX
        PARAMETER (NECUAMAX=(NMAPMAX+1)*NMAPMAX)
        INTEGER MXCOR
        PARAMETER (MXCOR=8) !numero maximo de intersecciones de dos pixels
!
        INTEGER READI,READILIM
        REAL READF
        REAL FPOLY
        CHARACTER*255 READC
!
        INTEGER J,J1,J2,I,I1,I2,K,K_,II,JJ
        INTEGER L2,J1_,J2_
        INTEGER IFINE,JFINE,IOPT
        INTEGER NCBUFF
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        INTEGER NMAP
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER IXMIN,IXMAX,IYMIN,IYMAX
        INTEGER NWIDTH,NMED
        INTEGER NEXTINFO
        INTEGER NF,NPEAKS
        INTEGER NDEG
        INTEGER NCOADD
        INTEGER NEWBUFF,NEWSCAN,NCOR
        INTEGER NSUBPIX,NSUBPIXMAX
        INTEGER NSPLFIT
        INTEGER NKNOTS
!!!        INTEGER I0BINS
        REAL AIJ(NECUAMAX),BIJ(NECUAMAX)
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL COEFFLOCAL(20),CHISQR
        REAL COEFF0(20),COEFF1(20),COEFF2(20),COEFF3(20)
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL XPEAK(NXYMAX)
        REAL YPEAK0(NXYMAX)
        REAL YPEAK1(NXYMAX),YPEAK2(NXYMAX),YPEAK3(NXYMAX)
        REAL EYPEAK1(NXYMAX),EYPEAK2(NXYMAX),EYPEAK3(NXYMAX)
        REAL YOPTIMAL(NXYMAX),XOPTIMAL(NXYMAX)
         REAL U,V,VINI,U1,U2,V1,V2
        REAL U1FINE,U2FINE,V1FINE,V2FINE
        REAL X1,Y1,X2,Y2,X3,Y3,X4,Y4
        REAL UMIN,UMAX,XMIN,XMAX,XMINMAX
        REAL YMIN,YMAX,DY,DX
        REAL XC,YC
        REAL XCINI,YCINI
        REAL X0,SIGMA,AMP,Y0
        REAL EX0,ESIGMA,EAMP
        REAL EEX0,EESIGMA,EEAMP,EEY0
        REAL FACTOR
        REAL XORIG,YORIG
        REAL XCOR1(2,4),XCOR2(2,4),XCOM(2,MXCOR)
        REAL FAREA,FAREATOT,FACTORFLUX
!!!        REAL SPLA(NXYMAX),SPLB(NXYMAX),SPLC(NXYMAX),SPLS(NXYMAX)
        REAL XKNOT(2),YRMSTOL,FRMS
        CHARACTER*1 CH,COPTIMAL
        CHARACTER*1 CREFINE,CMORE,CFIT
        CHARACTER*50 CDUMMY
        LOGICAL LFIT(NXYMAX)
        LOGICAL LOOP,IFU(NXMAX)
        LOGICAL LERR
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
        NWIDTH=11
!
        WRITE(*,101) 'Select spectrum point...'
        CALL PGSCI(5)
        CALL RPGBAND(0,0,0.,0.,XCINI,YCINI,CH)
        CALL PGSCI(1)
        IF(CH.EQ.'X')THEN
          RETURN
        END IF
!
        UMIN=REAL(SXMINGRID-SXMINEXTG)
        UMAX=REAL(SXMAXGRID+SXMAXEXTG)
        CALL FMAP(NMAP,AIJ,BIJ,XCINI,YCINI,U,VINI)
        DO J=1,NXYMAX
          U=UMIN+REAL(J-1)/REAL(NXYMAX-1)*(UMAX-UMIN)
          CALL FMAP(NMAP,AIJ_,BIJ_,U,VINI,XP(J),YP(J))
        END DO
        CALL PGSCI(5)
        CALL PGLINE(NXYMAX,XP,YP)
        CALL POLFIT(XP,YP,YP,NXYMAX,NMAP+1,0,COEFFLOCAL,CHISQR)
!
        CREFINE(1:1)=READC('Refine spectrum location','n','yn')
!------------------------------------------------------------------------------
        IF(CREFINE.EQ.'y')THEN
          WRITE(*,101) 'Select X-range...'
          CALL RPGBAND(6,0,0.,0.,XMIN,YC,CH)
          CALL RPGBAND(4,0,XMIN,0.,XMAX,YC,CH)
          CALL PGSCI(1)
          IF(XMIN.GT.XMAX)THEN
            XMINMAX=XMIN
            XMIN=XMAX
            XMAX=XMINMAX
          END IF
          DO J=1,NXYMAX
             XP(J)=XMIN+REAL(J-1)/REAL(NXYMAX-1)*(XMAX-XMIN)
             YP(J)=FPOLY(NMAP,COEFFLOCAL,XP(J))
          END DO
          CALL PGSCI(4)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
!
          CMORE(1:1)=READC('Plot individual fits (y/n)','n','yn')
          WRITE(*,101) '* Fit types:'
          WRITE(*,101) '  (1) GAUSSFIT'
          WRITE(*,101) '  (2) GAUSSFITA'
          WRITE(*,101) '  (3) GAUSCFIT'
          WRITE(*,101) '  (4) = (1) + (2)'
          WRITE(*,101) '  (5) = (2) + (3)'
          WRITE(*,101) '  (6) = (1) + (3)'
          WRITE(*,101) '  (7) = (1) + (2) + (3)'
          CFIT(1:1)=READC('Option','1','1234567')
!
          WRITE(CDUMMY,*) NWIDTH
          NWIDTH=READI('NWIDTH (must be odd)',CDUMMY)
!
          NMED=NWIDTH/2
          J1=NINT(XMIN)
          J2=NINT(XMAX)
          NEXTINFO=0
          DO J=J1,J2
            XC=REAL(J)
            YC=FPOLY(NMAP,COEFFLOCAL,XC)
            XPEAK(J-J1+1)=XC
            YPEAK0(J-J1+1)=YC
            I1=NINT(YC)-NMED
            I2=NINT(YC)+NMED
            IF(I1.LT.1)THEN
              I1=1
              I2=I1+NWIDTH-1
            ELSEIF(I2.GT.NAXIS(2,NCBUFF))THEN
              I2=NAXIS(2,NCBUFF)
              I1=I2-NWIDTH+1
            END IF
            DO I=I1,I2
              XF(I-I1+1)=REAL(I)
              YF(I-I1+1)=IMAGEN(J,I,NCBUFF)
            END DO
            NF=I2-I1+1
            IF(CMORE.EQ.'y')THEN
              CALL SUBPLOT(NF,1,NF,XF,YF,XF,YF,.TRUE.,.TRUE.,.FALSE.,.FALSE.,'y-axis','signal',' ',3,201,1.0)
            END IF
!
              IF(INDEX('1467',CFIT).NE.0)THEN
              CALL GAUSSFIT(NF,XF,YF,YF,X0,SIGMA,AMP,EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,1.E-6,0)
              YPEAK1(J-J1+1)=X0-YPEAK0(J-J1+1)
              EYPEAK1(J-J1+1)=EEX0
              IF(CMORE.EQ.'y')THEN
                DO I=1,NXYMAX
                  XP(I)=REAL(I1)+REAL(I-1)/REAL(NXYMAX-1)*REAL(I2-I1)
                  FACTOR=(XP(I)-X0)*(XP(I)-X0)/(2.*SIGMA*SIGMA)
                  IF(FACTOR.GT.60.)THEN
                    YP(I)=0.
                  ELSE
                    YP(I)=AMP*EXP(-FACTOR)
                  END IF
                END DO
                CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XP,YP,XP,YP,.FALSE.,.FALSE.,2,101,1.0)
                print*,'x0...: ',x0,ex0,eex0
                print*,'sigma: ',sigma,esigma,eesigma
                print*,'amp..: ',amp,eamp,eeamp
              END IF
            END IF
!
            IF(INDEX('2457',CFIT).NE.0)THEN
              CALL GAUSSFITA(NF,XF,YF,YF,X0,SIGMA,AMP,EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,1.E-6,0)
              YPEAK2(J-J1+1)=X0-YPEAK0(J-J1+1)
              EYPEAK2(J-J1+1)=EEX0
              IF(CMORE.EQ.'y')THEN
                DO I=1,NXYMAX
                  XP(I)=REAL(I1)+REAL(I-1)/REAL(NXYMAX-1)*REAL(I2-I1)
                  FACTOR=(XP(I)-X0)*(XP(I)-X0)/(2.*SIGMA*SIGMA)
                  IF(FACTOR.GT.60.)THEN
                    YP(I)=0.
                  ELSE
                    YP(I)=AMP*EXP(-FACTOR)
                  END IF
                END DO
                CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XP,YP,XP,YP,.FALSE.,.FALSE.,6,101,1.0)
                print*,'x0...: ',x0,ex0,eex0
                print*,'sigma: ',sigma,esigma,eesigma
                print*,'amp..: ',amp,eamp,eeamp
              END IF
            END IF
!
            IF(INDEX('3567',CFIT).NE.0)THEN
              CALL GAUSCFIT(NF,XF,YF,X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,1.E-6)
              YPEAK3(J-J1+1)=X0-YPEAK0(J-J1+1)
              EYPEAK3(J-J1+1)=EEX0
              IF(CMORE.EQ.'y')THEN
                DO I=1,NXYMAX
                  XP(I)=REAL(I1)+REAL(I-1)/REAL(NXYMAX-1)*REAL(I2-I1)
                  FACTOR=(XP(I)-X0)*(XP(I)-X0)/(2.*SIGMA*SIGMA)
                  IF(FACTOR.GT.60.)THEN
                    YP(I)=Y0
                  ELSE
                    YP(I)=Y0+AMP*EXP(-FACTOR)
                  END IF
                END DO
                CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XP,YP,XP,YP,.FALSE.,.FALSE.,5,101,1.0)
                print*,'x0...: ',x0,eex0
                print*,'sigma: ',sigma,eesigma
                print*,'amp..: ',amp,eeamp
                print*,'y0...: ',y0,eey0
              END IF
            END IF
            IF(CMORE.EQ.'y')THEN
              CMORE(1:1)=READC('More plots (y/n)','y','yn')
            END IF
            CALL SHOWPERC(J1,J2,1,J,NEXTINFO)
          END DO
          NPEAKS=J2-J1+1
          CMORE='y'
          DO WHILE(CMORE.EQ.'y')
            IF(INDEX('1467',CFIT).NE.0)THEN
              CALL SUBPLOT(NPEAKS,1,NPEAKS,XPEAK,YPEAK1,XPEAK,EYPEAK1,.TRUE.,.TRUE.,.FALSE.,.TRUE.,'X axis','\\gDY',' ',3,1,1.0)
            END IF
            IF(INDEX('2457',CFIT).NE.0)THEN
              IF(INDEX('47',CFIT).NE.0)THEN
                CALL SUBPLOTBIS(NPEAKS,1,NPEAKS,XPEAK,YPEAK2,XPEAK,EYPEAK2,.FALSE.,.TRUE.,4,1,1.0)
              ELSE
                CALL SUBPLOT(NPEAKS,1,NPEAKS,XPEAK,YPEAK2,XPEAK,EYPEAK2,.TRUE.,.TRUE.,.FALSE.,.TRUE.,'X axis','\\gDY',' ',3,1,1.0)
              END IF
            END IF
            IF(INDEX('3567',CFIT).NE.0)THEN
              IF(CFIT.EQ.'3')THEN
                CALL SUBPLOT(NPEAKS,1,NPEAKS,XPEAK,YPEAK3,XPEAK,EYPEAK3,.TRUE.,.TRUE.,.FALSE.,.TRUE.,'X axis','\\gDY',' ',3,1,1.0)
              ELSE
                CALL SUBPLOTBIS(NPEAKS,1,NPEAKS,XPEAK,YPEAK3,XPEAK,EYPEAK3,.FALSE.,.TRUE.,5,1,1.0)
              END IF
             END IF
            NDEG=READILIM('Polynomial degree','0',0,MIN0(19,NPEAKS-1))
            IF(INDEX('1467',CFIT).NE.0)THEN
              CALL POLFITSIG(NPEAKS,XPEAK,YPEAK1,3.0,NDEG,COEFF1,LFIT)
              DO J=1,NXYMAX
                XP(J)=XMIN+REAL(J-1)/REAL(NXYMAX-1)*(XMAX-XMIN)
                YP(J)=FPOLY(NDEG,COEFF1,XP(J))
              END DO
              CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XP,YP,XP,YP,.FALSE.,.FALSE.,2,101,1.0)
            END IF
            IF(INDEX('2457',CFIT).NE.0)THEN
              CALL POLFITSIG(NPEAKS,XPEAK,YPEAK2,3.0,NDEG,COEFF2,LFIT)
              DO J=1,NXYMAX
                XP(J)=XMIN+REAL(J-1)/REAL(NXYMAX-1)*(XMAX-XMIN)
                YP(J)=FPOLY(NDEG,COEFF2,XP(J))
              END DO
              CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XP,YP,XP,YP,.FALSE.,.FALSE.,6,101,1.0)
            END IF
            IF(INDEX('3567',CFIT).NE.0)THEN
              CALL POLFITSIG(NPEAKS,XPEAK,YPEAK3,3.0,NDEG,COEFF3,LFIT)
              DO J=1,NXYMAX
                XP(J)=XMIN+REAL(J-1)/REAL(NXYMAX-1)*(XMAX-XMIN)
                YP(J)=FPOLY(NDEG,COEFF3,XP(J))
              END DO
              CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XP,YP,XP,YP,.FALSE.,.FALSE.,7,101,1.0)
            END IF
            CMORE(1:1)=READC('Another fit (y/n)','n','yn')
          END DO
!
          IF(CFIT.EQ.'4')THEN
            CFIT(1:1)=READC('Fit (0=None/1/2)','@','012')
          ELSEIF(CFIT.EQ.'5')THEN
            CFIT(1:1)=READC('Fit (0=None/2/3)','@','023')
          ELSEIF(CFIT.EQ.'6')THEN
            CFIT(1:1)=READC('Fit (0=None/1/3)','@','013')
          ELSEIF(CFIT.EQ.'7')THEN
            CFIT(1:1)=READC('Fit (0=None/1/2/3)','@','0123')
          END IF
!
          IF(CFIT.EQ.'0')THEN
            NDEG=0
            COEFF0(1)=0.
          ELSEIF(CFIT.EQ.'1')THEN
            DO K=1,NDEG+1
              COEFF0(K)=COEFF1(K)
            END DO
          ELSEIF(CFIT.EQ.'2')THEN
            DO K=1,NDEG+1
              COEFF0(K)=COEFF2(K)
            END DO
          ELSEIF(CFIT.EQ.'3')THEN
            DO K=1,NDEG+1
              COEFF0(K)=COEFF3(K)
            END DO
          END IF
        ELSE
          NDEG=0
          COEFF0(1)=0.
        END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! NDEG,COEFF0 es el polinomio correctivo que aplicamos a NMAP,COEFFLOCAL para
! determinar el centro del perfil espacial.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Dibujamos el nuevo centro del objeto
        DO J=1,NXYMAX
          U=UMIN+REAL(J-1)/REAL(NXYMAX-1)*(UMAX-UMIN)
          CALL FMAP(NMAP,AIJ_,BIJ_,U,VINI,XP(J),YP(J))
          YP(J)=YP(J)+FPOLY(NDEG,COEFF0,XP(J))
        END DO
        CALL PGSCI(2)
        CALL PGLINE(NXYMAX,XP,YP)
        CALL PGSCI(1)
! NCOADD se refiere a pixels en el sistema distorsionado, no en la matriz
! original de datos.
        NCOADD=READI('Number of spectra at each side of maximum','@')
! Dibujamos la region a sumar
        DO J=1,NXYMAX
          U=UMIN+REAL(J-1)/REAL(NXYMAX-1)*(UMAX-UMIN)
          CALL FMAP(NMAP,AIJ_,BIJ_,U,VINI-NCOADD,XP(J),YP(J))
          YP(J)=YP(J)+FPOLY(NDEG,COEFF0,XP(J))
        END DO
        CALL PGSCI(7)
        CALL PGLINE(NXYMAX,XP,YP)
        DO J=1,NXYMAX
          U=UMIN+REAL(J-1)/REAL(NXYMAX-1)*(UMAX-UMIN)
          CALL FMAP(NMAP,AIJ_,BIJ_,U,VINI+NCOADD,XP(J),YP(J))
          YP(J)=YP(J)+FPOLY(NDEG,COEFF0,XP(J))
        END DO
        CALL PGLINE(NXYMAX,XP,YP)
        CALL PGSCI(1)
        NEWBUFF=READILIM('New buffer# to store spectrum','@',1,NMAXBUFF)
        NEWSCAN=READILIM('New scan# to store spectrum','@',1,NYMAX)
        COPTIMAL(1:1)=READC('Optimal extraction (y/n)','n','yn')
        LERR=(NCBUFF.LE.NMAXBUFF/2).AND.(NEWBUFF.LE.NMAXBUFF/2)
        J1=SXMINGRID-SXMINEXTG
        J2=SXMAXGRID+SXMAXEXTG
        IF(NAXIS(1,NEWBUFF).LT.J2-J1) NAXIS(1,NEWBUFF)=J2-J1
        IF(NAXIS(2,NEWBUFF).LT.NEWSCAN) NAXIS(2,NEWBUFF)=NEWSCAN
        IF(LERR)THEN
          IF(NAXIS(1,NEWBUFF+NMAXBUFF/2).LT.J2-J1) NAXIS(1,NEWBUFF+NMAXBUFF/2)=J2-J1
          IF(NAXIS(2,NEWBUFF+NMAXBUFF/2).LT.NEWSCAN) NAXIS(2,NEWBUFF+NMAXBUFF/2)=NEWSCAN
        END IF
!------------------------------------------------------------------------------
! Extraccion no optimizada
!------------------------------------------------------------------------------
        IF(COPTIMAL.EQ.'n')THEN
          IF(NAXIS(1,NEWBUFF).LT.J2-J1) NAXIS(1,NEWBUFF)=J2-J1
          IF(NAXIS(2,NEWBUFF).LT.NEWSCAN) NAXIS(2,NEWBUFF)=NEWSCAN
          NEXTINFO=0
          DO J=J1,J2-1
            IMAGEN(J-J1+1,NEWSCAN,NEWBUFF)=0.
            IF(LERR) IMAGEN(J-J1+1,NEWSCAN,NEWBUFF+NMAXBUFF/2)=0.
            U1=REAL(J) !borde izquierdo
            U2=U1+1.0  !borde derecho
            CALL FMAP(NMAP,AIJ_,BIJ_,(U1+U2)/2.,VINI,XC,YC)
            YC=YC+FPOLY(NDEG,COEFF0,XC)
            CALL FMAP(NMAP,AIJ,BIJ,XC,YC,U,V)
            DO I=-NCOADD,NCOADD
              V1=V-0.5+I
              V2=V1+1.0
              CALL FMAP(NMAP,AIJ_,BIJ_,U1,V1,X1,Y1)
              CALL FMAP(NMAP,AIJ_,BIJ_,U2,V1,X2,Y2)
              CALL FMAP(NMAP,AIJ_,BIJ_,U2,V2,X3,Y3)
              CALL FMAP(NMAP,AIJ_,BIJ_,U1,V2,X4,Y4)
!!!              IF(I.EQ.0)THEN
!!!                CALL PPOLY4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,2)
!!!              ELSE
!!!                CALL PPOLY4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,3)
!!!              END IF
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
              DO II=IYMIN,IYMAX
                DO JJ=IXMIN,IXMAX
                  IF((II.GE.1).AND.(II.LE.NAXIS(2,NCBUFF)).AND.(JJ.GE.1).AND.(JJ.LE.NAXIS(1,NCBUFF)))THEN
                    X1=REAL(JJ)-0.5
                    X2=X1+1.
                    Y1=REAL(II)-0.5
                    Y2=Y1+1.
!!!                 CALL PPOLY4(X1,Y1,X2,Y1,X2,Y2,X1,Y2,6)
                    XCOR1(1,1)=X1-XORIG
                    XCOR1(2,1)=Y1-YORIG
                    XCOR1(1,2)=X2-XORIG
                    XCOR1(2,2)=XCOR1(2,1)
                    XCOR1(1,3)=XCOR1(1,2)
                    XCOR1(2,3)=Y2-YORIG
                    XCOR1(1,4)=XCOR1(1,1)
                    XCOR1(2,4)=XCOR1(2,3)
!!!                 CALL PPOLYN(4,XCOR1,FAREATOT)
                    FAREATOT=1.0 !lo sabemos, por lo que no hacemos PPOLYN
                    CALL SRCPAT(4,XCOR1,4,XCOR2,NCOR,XCOM)
                    IF(NCOR.GT.0)THEN
                      CALL PPOLYN(NCOR,XCOM,FAREA)
                      IMAGEN(J-J1+1,NEWSCAN,NEWBUFF)=IMAGEN(J-J1+1,NEWSCAN,NEWBUFF)+IMAGEN(JJ,II,NCBUFF)*FAREA/FAREATOT
                      IF(LERR)THEN
                        IMAGEN(J-J1+1,NEWSCAN,NEWBUFF+NMAXBUFF/2)=IMAGEN(J-J1+1,NEWSCAN,NEWBUFF+NMAXBUFF/2)+ &
                         IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)*IMAGEN(JJ,II,NCBUFF+NMAXBUFF/2)*FAREA/FAREATOT
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END DO
            IF(LERR)THEN
              IMAGEN(J-J1+1,NEWSCAN,NEWBUFF+NMAXBUFF/2)=SQRT(IMAGEN(J-J1+1,NEWSCAN,NEWBUFF+NMAXBUFF/2))
            END IF
            CALL SHOWPERC(J1,J2-1,1,J,NEXTINFO)
          END DO
          RETURN
        END IF
!------------------------------------------------------------------------------
! Extraccion optimizada
        WRITE(*,101) '****************************************'
        WRITE(*,101) 'Note that this method has not been fully'
        WRITE(*,101) 'implemented and ERRORS have not been'
        WRITE(*,101) 'included so far.'
        WRITE(*,101) '****************************************'
        WRITE(*,100) 'Press <CR> to continue...'
        READ(*,*)
!------------------------------------------------------------------------------
! Definimos el recorrido espectral a utilizar (usamos coordenadas en el sistema
! de referencia del espectro distorsionado; es lo que en el programa llamo
! U range).
        DO J=1,NXMAX
          IFU(J)=.FALSE.
        END DO
        WRITE(CDUMMY,*) J2 !OJO: es J2 y no J2-1 porque no existe IFU(0)
        CALL RMBLANK(CDUMMY,CDUMMY,L2)
        LOOP=.TRUE.
        K=0
        DO WHILE(LOOP)
          IF(K.EQ.0)THEN
            CALL READ2I('U range (0,0=exit)','1,'//CDUMMY(1:L2),J1_,J2_)
          ELSE
            CALL READ2I('U range (0,0=exit)','0,0',J1_,J2_)
          END IF
          IF((J1_.EQ.0).AND.(J2_.EQ.0))THEN
            K=0
            DO J=J1,J2-1
              IF(IFU(J+1)) K=K+1 !OJO: es IFU(J+1) porque no existe IFU(0)
            END DO
            IF(K.GT.0)THEN
              LOOP=.FALSE.
            ELSE
              WRITE(*,101) 'ERROR: No U range selected. Try again.'
            END IF
          ELSE
            IF(J1_.LT.J1+1) J1_=J1+1
            IF(J2_.GT.J2) J2_=J2
            DO J=J1_,J2_
              IFU(J)=.TRUE.
            END DO
            DO J=1,NXYMAX
              U=J1_-1+REAL(J-1)/REAL(NXYMAX-1)*(J2_-J1_+1)
              CALL FMAP(NMAP,AIJ_,BIJ_,U,VINI-NCOADD,XP(J),YP(J))
              YP(J)=YP(J)+FPOLY(NDEG,COEFF0,XP(J))
            END DO
            CALL PGSCI(6)
            CALL PGLINE(NXYMAX,XP,YP)
            DO J=1,NXYMAX
              U=J1_-1+REAL(J-1)/REAL(NXYMAX-1)*(J2_-J1_+1)
              CALL FMAP(NMAP,AIJ_,BIJ_,U,VINI+NCOADD,XP(J),YP(J))
              YP(J)=YP(J)+FPOLY(NDEG,COEFF0,XP(J))
            END DO
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
            K=0
            DO J=J1,J2-1
              IF(IFU(J+1)) K=K+1 !OJO: es IFU(J+1) porque no existe IFU(0)
            END DO
          END IF
        END DO
!..............................................................................
! Definimos el subpixelado a emplear. Este factor es tenido en cuenta por
! la variable FACTORFLUX que garantiza que el corte promedio tenga unidades
! similares a un corte promedio de la imagen.
        NSUBPIXMAX=NXYMAX/(2*NCOADD+1)        !OJO: division con truncamiento
        NSUBPIX=READILIM('No. of subpixels','1',1,NSUBPIXMAX)
        NSPLFIT=(2*NCOADD+1)*NSUBPIX
        FACTORFLUX=REAL(NSUBPIX)/REAL(K)
!..............................................................................
! Las variables XOPTIMAL,YOPTIMAL van a almacenar el corte promedio en la
! direccion espacial con el subpixelado elegido.
        DO IOPT=1,NSPLFIT
          XOPTIMAL(IOPT)=-REAL(NCOADD)-0.5+(REAL(IOPT)-0.5)/REAL(NSUBPIX)
          YOPTIMAL(IOPT)=0.
        END DO
!..............................................................................
! Calculamos el corte espacial promedio.
        NEXTINFO=0
        K_=0
        DO J=J1,J2-1
          IF(IFU(J+1))THEN !OJO: es IFU(J+1) porque no existe IFU(0)
            K_=K_+1
            U1=REAL(J) !borde izquierdo pixel original
            DO JFINE=1,NSUBPIX
              U1FINE=U1+REAL(JFINE-1)/REAL(NSUBPIX) !borde izquierdo subpixel
              U2FINE=U1FINE+1.0/REAL(NSUBPIX)       !borde derecho subpixel
              CALL FMAP(NMAP,AIJ_,BIJ_,(U1FINE+U2FINE)/2.,VINI,XC,YC)
              YC=YC+FPOLY(NDEG,COEFF0,XC)
              CALL FMAP(NMAP,AIJ,BIJ,XC,YC,U,V)
              IOPT=0
              DO I=-NCOADD,NCOADD
                V1=V-0.5+I !borde inferior pixel original
                DO IFINE=1,NSUBPIX
                  IOPT=IOPT+1
                  V1FINE=V1+REAL(IFINE-1)/REAL(NSUBPIX)!borde inferior subpixel
                  V2FINE=V1FINE+1./REAL(NSUBPIX)       !borde superior subpixel
                  CALL FMAP(NMAP,AIJ_,BIJ_,U1FINE,V1FINE,X1,Y1)
                  CALL FMAP(NMAP,AIJ_,BIJ_,U2FINE,V1FINE,X2,Y2)
                  CALL FMAP(NMAP,AIJ_,BIJ_,U2FINE,V2FINE,X3,Y3)
                  CALL FMAP(NMAP,AIJ_,BIJ_,U1FINE,V2FINE,X4,Y4)
!!!                  IF(I.EQ.0)THEN
!!!                    CALL PPOLY4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,2)
!!!                  ELSE
!!!                    CALL PPOLY4(X1,Y1,X2,Y2,X3,Y3,X4,Y4,3)
!!!                  END IF
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
                  DO II=IYMIN,IYMAX
                    DO JJ=IXMIN,IXMAX
                      IF((II.GE.1).AND.(II.LE.NAXIS(2,NCBUFF)).AND.(JJ.GE.1).AND.(JJ.LE.NAXIS(1,NCBUFF)))THEN
                        X1=REAL(JJ)-0.5
                        X2=X1+1.
                        Y1=REAL(II)-0.5
                        Y2=Y1+1.
!!!                        CALL PPOLY4(X1,Y1,X2,Y1,X2,Y2,X1,Y2,6)
                        XCOR1(1,1)=X1-XORIG
                        XCOR1(2,1)=Y1-YORIG
                        XCOR1(1,2)=X2-XORIG
                        XCOR1(2,2)=XCOR1(2,1)
                        XCOR1(1,3)=XCOR1(1,2)
                        XCOR1(2,3)=Y2-YORIG
                        XCOR1(1,4)=XCOR1(1,1)
                        XCOR1(2,4)=XCOR1(2,3)
!!!                        CALL PPOLYN(4,XCOR1,FAREATOT)
                        FAREATOT=1.0 !lo sabemos, por lo que no hacemos PPOLYN
                        CALL SRCPAT(4,XCOR1,4,XCOR2,NCOR,XCOM)
                        IF(NCOR.GT.0)THEN
                          CALL PPOLYN(NCOR,XCOM,FAREA)
                          YOPTIMAL(IOPT)=YOPTIMAL(IOPT)+IMAGEN(JJ,II,NCBUFF)*FAREA/FAREATOT*FACTORFLUX
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END IF
          CALL SHOWPERC(1,K,1,K_,NEXTINFO)
        END DO
        LOOP=.TRUE.
        CALL FINDMML(NSPLFIT,1,NSPLFIT,XOPTIMAL,XMIN,XMAX)
        CALL FINDMML(NSPLFIT,1,NSPLFIT,YOPTIMAL,YMIN,YMAX)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/20.
        XMAX=XMAX+DX/20.
        DY=YMAX-YMIN
        YMIN=YMIN-DY/20.
        YMAX=YMAX+DY/20.
        DO WHILE(LOOP)
          CALL RPGERASW(0.4,1.0,0.0,0.8,0)
          CALL PGSVP(0.5,0.95,0.10,0.75)
          CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
          CALL PGLABEL('spatial direction','No. of counts',' ')
          CALL PGSCI(3)
          CALL PGBIN(NSPLFIT,XOPTIMAL,YOPTIMAL,.TRUE.)
          CALL PGSCI(1)
          NKNOTS=READILIM('No. of knots','8',3,20)
          DO I=1,NKNOTS
            XKNOT(I)=XOPTIMAL(1)+REAL(I-1)/REAL(NKNOTS-1)*(XOPTIMAL(NSPLFIT)-XOPTIMAL(1))
          END DO
          YRMSTOL=READF('YRMSTOL for downhill','1.E-6')
          CALL PGSCI(2)
          CALL SPLFIT(NSPLFIT,XOPTIMAL,YOPTIMAL,NKNOTS,XKNOT,YRMSTOL,NXYMAX,XP,YP,XOPTIMAL(1),XOPTIMAL(NSPLFIT),FRMS,.TRUE.)
          CALL PGSCI(1)
          CMORE(1:1)=READC('Another fit (y/n)','n','yn')
          LOOP=(CMORE.EQ.'y')
        END DO
!..............................................................................
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
