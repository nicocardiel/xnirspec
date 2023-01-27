!------------------------------------------------------------------------------
! Version 27-July-1998                                           File: splfit.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!omment
!
! SUBROUTINE SPLFIT(N,X,Y,ND,XD,YRMSTOL,NOUT,XOUT,YOUT,XMIN,XMAX,SIGMA,LPLOTS)
!
! Input: N,X,Y,ND,XD,YRMSTOL,NOUT,XMIN,XMAX,SIGMA,LPLOTS
! Output: XOUT,YOUT
!
! Least-squares fit to splines, using ND knots located at XD(). 
! Input data are X(N), Y(N). XOUT(NOUT), YOUT(NOUT) are the output values which
! are computed in the range from XMIN to XMAX. The knot location determines the
! X(),Y() range employed in the fit (which is performed in the interval from
! XD(1) to XD(ND))
!
! INTEGER N -> initial number of points in input data
! REAL    X(N) -> sorted input data
! REAL    Y(N) -> input data
! INTEGER ND -> number of knots
! REAL    XD(ND) -> X location of each knot
! REAL    YRMSTOL ->  stopping criterion for DOWNHILL
! INTEGER NOUT -> number of points in output
! REAL    XOUT(NOUT) -> output data
! REAL    YOUT(NOUT) -> output data
! REAL    XMIN -> = XOUT(1)
! REAL    XMAX -> = XOUT(NOUT)
! REAL    SIGMA -> sigma of the fit
! LOGICAL LPLOTS -> if .TRUE. some plots are performed
!
!omment
!------------------------------------------------------------------------------
! Esta subrutina requiere POLFIT y DOWNHILL
        SUBROUTINE SPLFIT(N,X,Y,ND,XD,YRMSTOL,NOUT,XOUT,YOUT,XMIN,XMAX,SIGMA,LPLOTS)
        IMPLICIT NONE
        INCLUDE 'largest.inc'
! NDMAX- numero maximo de nodos
        INTEGER NDMAX                         !cambiar tambien en funcion FUNK
        PARAMETER (NDMAX=20)
!
        INTEGER READI,READILIM
        REAL READF
        CHARACTER*255 READC
!
        INTEGER N
        REAL X(N),Y(N)
        INTEGER ND
        REAL XD(ND)
        INTEGER NOUT
        REAL XOUT(NOUT),YOUT(NOUT)
        REAL XMIN,XMAX,SIGMA
        LOGICAL LPLOTS
!
        INTEGER I,J,K,H,L
        INTEGER NF
        INTEGER NEVAL
        INTEGER NDD,NREF
        INTEGER NITER,NITERT
        INTEGER NCOLOR
        INTEGER NSEED,NRANND(NDMAX)
        INTEGER I0SPL
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL YD(NDMAX),MEANF,RMSF
        REAL ASPL(NDMAX),BSPL(NDMAX),CSPL(NDMAX),SSPL(NDMAX)
        REAL XDD(NDMAX),XX(NDMAX),DXX(NDMAX),XX0(NDMAX),DXX0(NDMAX)
        REAL YRMSTOL
        REAL YINI,YFIN
        EXTERNAL FUNKSPLFIT,FUNKSPLFIT1,FUNKSPLFIT2,FUNKSPLFIT3
        REAL FUNKSPLFIT,FUNKSPLFIT1,FUNKSPLFIT2,FUNKSPLFIT3
        REAL A(4),CHISQR,FDUMMY
        REAL YMINP,YMAXP,XMINP,XMAXP
        CHARACTER*1 CREF
        CHARACTER*20 CDUMMY
!
        COMMON/BLKNSEED/NSEED
        COMMON/BLKFUNK1/NF
        COMMON/BLKFUNK2/XF,YF
        COMMON/BLKFUNK3/NDD
        COMMON/BLKFUNK4/XDD
        COMMON/BLKFUNK5/YD
        COMMON/BLKFUNK6/NREF
!------------------------------------------------------------------------------
        NITER=0
        NITERT=0
        YINI=0 !evita WARNING de compilacion
        YFIN=0 !evita WARNING de compilacion
        CREF='X' !evita WARNING de compilacion
        if(lplots)then
          call pgqwin(xminp,xmaxp,yminp,ymaxp)
          call pgqci(ncolor)
        end if
        IF(ND.GT.NDMAX)THEN
          WRITE(*,*)
          WRITE(*,101)'ERROR in SPLFIT: '
          WRITE(*,110)'>>> No. of Knots: ',ND
          WRITE(*,110)'>>> Maximum No. of Knots: ',NDMAX
          WRITE(*,*)
          GOTO 900
        END IF
        DO I=1,ND
          YD(I)=0.
        END DO
        DO I=1,NOUT
          XOUT(I)=XMIN+(XMAX-XMIN)*REAL(I-1)/REAL(NOUT-1)
        END DO
! ajuste inicial a polinomio de grado 3
        J=1
        DO WHILE(X(J).LT.XD(1))
          J=J+1
        END DO
        I=2
        DO WHILE(I.LE.ND)
          K=1
          DO WHILE((X(J).LE.XD(I)).AND.(J.LE.N))
            XF(K)=X(J)
            YF(K)=Y(J)
            K=K+1
            IF(K.GT.NXYMAX)THEN
              WRITE(*,*)
              WRITE(*,101)'ERROR in SPLFIT:'
              WRITE(*,101)'>>> No. of points to fit > NXYMAX'
              WRITE(*,*)
              GOTO 900
            END IF
            J=J+1
          END DO
          IF(K.LT.4)THEN            !tomamos el valor medio de los puntos si no
            FDUMMY=0.               !hay suficientes para el ajuste polinomico
            DO H=NINT(XD(I-1)),NINT(XD(I))
              FDUMMY=FDUMMY+YF(H)
            END DO
            FDUMMY=FDUMMY/REAL(NINT(XD(I))-NINT(XD(I-1))+1)
            YD(I-1)=YD(I-1)+FDUMMY
            YD(I)=YD(I)+FDUMMY
          ELSE
            NF=K-1
            CALL POLFIT(XF,YF,YF,NF,4,0,A,CHISQR)
            YD(I-1)=YD(I-1)+A(1)+A(2)*XD(I-1)+A(3)*XD(I-1)*XD(I-1)+A(4)*XD(I-1)*XD(I-1)*XD(I-1)
            YD(I)=YD(I)+A(1)+A(2)*XD(I)+A(3)*XD(I)*XD(I)+A(4)*XD(I)*XD(I)*XD(I)
! dibujamos polinomios iniciales de grado 3
!!!            if(lplots)then
!!!              do k=1,100
!!!                xf(k)=xd(i-1)+(xd(i)-xd(i-1))*real(k-1)/99.
!!!                yf(k)=a(1)+a(2)*xf(k)+a(3)*xf(k)**2+a(4)*xf(k)**3
!!!              end do
!!!              call pgsci(i+1)
!!!              call pgline(100,xf,yf)
!!!            end if
          END IF
          I=I+1
          IF(X(J-1).EQ.XD(I-1)) J=J-1
        END DO
        DO I=2,ND-1
          YD(I)=YD(I)/2.
        END DO
!
        J=1
        DO WHILE(X(J).LT.XD(1))
          J=J+1
        END DO
        I=1
        DO WHILE((X(J).LE.XD(ND)).AND.(J.LE.N))
          XF(I)=X(J)
          YF(I)=Y(J)
          I=I+1
          IF(I.GT.NXYMAX)THEN
            WRITE(*,*)
            WRITE(*,101)'ERROR in SPLFIT:'
            WRITE(*,101)'Please, redim NXYMAX'
            WRITE(*,*)
            GOTO 900
          END IF
          J=J+1
        END DO
        NF=I-1
!
        NDD=ND
        DO I=1,ND
          XDD(I)=XD(I)
        END DO
! calculamos media y rms de YD
        MEANF=0.
        DO I=1,ND
          MEANF=MEANF+YD(I)
        END DO
        MEANF=MEANF/REAL(ND)
        RMSF=0.
        DO I=1,ND
          RMSF=RMSF+(MEANF-YD(I))*(MEANF-YD(I))
        END DO
        RMSF=SQRT(RMSF/REAL(ND-1))
! valores iniciales para DOWNHILL
        DO I=1,ND
          XX0(I)=YD(I)
          IF(RMSF.GT.0.0)THEN
            DXX0(I)=RMSF/3.
          ELSEIF(YD(I).GT.0.0)THEN
            DXX0(I)=0.05*YD(I)
          ELSE
            DXX0(I)=1.0                        !que remedio; no tenemos ni idea
          END IF
        END DO
! llamamos a DOWNHILL
        WRITE(*,100)'Running DOWNHILL...'
        CALL DOWNHILL(ND,XX0,DXX0,FUNKSPLFIT,1.0,0.5,2.0,YRMSTOL,XX,DXX,NEVAL,500)
        WRITE(*,110)'      no. of function evaluations: ',NEVAL
        NSEED=-1       !semilla inicial para el generador de numeros aleatorios
        DO J=1,ND
          YD(J)=XX(J)
        END DO
        SIGMA=SQRT(FUNKSPLFIT(YD))
20      SSPL(1)=0.
        SSPL(ND)=0.
        CALL CUBSPL(XD,YD,ND,4,SSPL,ASPL,BSPL,CSPL)                    !IMODE=1
        I0SPL=1                  !la primera vez busca en el inicio de la tabla
        DO K=1,NOUT
          CALL CUBSPLX(XD,YD,ASPL,BSPL,CSPL,ND,I0SPL,XOUT(K),YOUT(K))
        END DO
! si estamos iterando seguimos con las iteraciones
        IF((NITERT.NE.0).AND.(NITER.LT.NITERT)) GOTO 24
!
! dibujamos ajuste final
        if(lplots)then
            if(ncolor.gt.1) call pgsci(1)
            call pgline(nout,xout,yout)
            do i=1,nd
              call pgpoint(1,xd(i),yd(i),17)
              write(cdummy,*)i
              call rmblank(cdummy,cdummy,k)
              call pgptext(xd(i),yd(i)+(ymaxp-yminp)/25.,0.,0.5,cdummy(1:k))
            end do
        end if
        write(*,100)'sigma of the fit: '
        write(*,*)sigma
!
! si el numero de Knots es solo 2 (los extremos) no se permite refinar el
! ajuste
        IF(ND.EQ.2) RETURN
! Si se quiere refinamos el ajuste
        CREF(1:1)=READC('Refine the fit (y/n)','n','yn')
        IF(CREF.EQ.'n')THEN
          if(lplots)then       !pintamos ajuste final en grueso para distinguir
              if(ncolor.gt.1) call pgsci(1)
              call pgslw(3)
              call pgline(nout,xout,yout)
              do i=1,nd
                call pgpoint(1,xd(i),yd(i),17)
                write(cdummy,*)i
                call rmblank(cdummy,cdummy,k)
                call pgptext(xd(i),yd(i)+(ymaxp-yminp)/25.,0.,0.5,cdummy(1:k))
              end do
              call pgslw(1)
          end if
          RETURN
        END IF
        WRITE(*,101) '(1) Refine X and Y position-> 1 Knot'
        WRITE(*,101) '(2) Refine X-position  -----> 1 Knot'
        WRITE(*,101) '(3) Refine Y-position  -----> 1 Knot'
        WRITE(*,100) '(W) Refine X and Y position-> all Knots'
        WRITE(*,101) ' (except first and last)'
        WRITE(*,101) '(0) Exit'
        CREF(1:1)=READC('Option','0','0123Ww')
        IF(CREF.EQ.'w')CREF='W'
        IF(CREF.EQ.'0')THEN
          if(lplots)then       !pintamos ajuste final en grueso para distinguir
              if(ncolor.gt.1) call pgsci(1)
              call pgslw(3)
              call pgline(nout,xout,yout)
              do i=1,nd
                call pgpoint(1,xd(i),yd(i),17)
                write(cdummy,*)i
                call rmblank(cdummy,cdummy,k)
                call pgptext(xd(i),yd(i)+(ymaxp-yminp)/25.,0.,0.5,cdummy(1:k))
              end do
              call pgslw(1)
          end if
          RETURN
        END IF
        NITER=0
        IF(CREF.NE.'W')THEN
          NREF=READILIM('Knot number to be refined','@',2,ND-1)
          NITERT=0
        ELSE
          YINI=READF('Function value in first knot','0.0')
          YFIN=READF('Function value in last  knot','0.0')
          NITERT=READI('How many iterations','1')
          NITER=0
        END IF
! pintamos con otro color la curva que va a quedar desfasada
        if(lplots)then
            if(ncolor.gt.1)then
              call pgsci(ncolor)
              call pgline(nout,xout,yout)
              do i=1,nd
                call pgpoint(1,xd(i),yd(i),17)
                write(cdummy,*)i
                call rmblank(cdummy,cdummy,k)
                call pgptext(xd(i),yd(i)+(ymaxp-yminp)/25.,0.,0.5,cdummy(1:k))
              end do
              call pgsci(1)
            end if
        end if
!
        WRITE(CDUMMY,*)YRMSTOL
        YRMSTOL=READF('YRMSTOL for downhill',CDUMMY)
! -> REFINAMOS x e y ----------------------------------------------------------
24      IF(CREF.EQ.'1')THEN
          WRITE(*,100)'Valor inicial  en X,Y: '
          WRITE(*,*) XD(NREF),YD(NREF)
          XX0(1)=XD(NREF)                      !valores iniciales para DOWNHILL
          IF(XD(NREF-1).NE.XD(NREF))THEN
            DXX0(1)=(XD(NREF)-XD(NREF-1))*0.05
          ELSEIF(XD(NREF+1).NE.XD(NREF))THEN
            DXX0(1)=(XD(NREF)-XD(NREF+1))*0.05
          ELSE
            DXX0(1)=1. !que remedio
          END IF
          XX0(2)=YD(NREF)
          IF(YD(NREF).NE.0.0)THEN
            DXX0(2)=YD(NREF)*0.05
          ELSEIF(YD(NREF-1).NE.YD(NREF))THEN
            DXX0(2)=(YD(NREF)-YD(NREF-1))*0.2
          ELSEIF(YD(NREF+1).NE.YD(NREF))THEN
            DXX0(2)=(YD(NREF)-YD(NREF+1))*0.2
          ELSE
            DXX0(2)=1. !que remedio
          END IF
          WRITE(*,100)'Running DOWNHILL...'
          CALL DOWNHILL(2,XX0,DXX0,FUNKSPLFIT3,1.0,0.5,2.0,YRMSTOL,XX,DXX,NEVAL,500)
          WRITE(*,110)'      no. of function evaluations: ',NEVAL
          XD(NREF)=XX(1)
          YD(NREF)=XX(2)
          WRITE(*,100)'Valor refinado en X,Y: '
          WRITE(*,*) XD(NREF),YD(NREF)
          SIGMA=SQRT(FUNKSPLFIT3(XX))
          DO I=1,ND            !actualizamos XDD para futuras llamadas a FUNK's
            XDD(I)=XD(I)
          END DO
! -> REFINAMOS x --------------------------------------------------------------
        ELSEIF(CREF.EQ.'2')THEN
          WRITE(*,100)'Valor inicial  en X: '
          WRITE(*,*) XD(NREF)
          XX0(1)=XD(NREF)                      !valores iniciales para DOWNHILL
          IF(XD(NREF-1).NE.XD(NREF))THEN
            DXX0(1)=(XD(NREF)-XD(NREF-1))*0.05
          ELSEIF(XD(NREF+1).NE.XD(NREF))THEN
            DXX0(1)=(XD(NREF)-XD(NREF+1))*0.05
          ELSE
            DXX0(1)=1. !que remedio
          END IF
          WRITE(*,100)'Running DOWNHILL...'
          CALL DOWNHILL(1,XX0,DXX0,FUNKSPLFIT1,1.0,0.5,2.0,YRMSTOL,XX,DXX,NEVAL,500)
          WRITE(*,110)'      no. of function evaluations: ',NEVAL
          XD(NREF)=XX(1)
          SIGMA=SQRT(FUNKSPLFIT1(XX))
          DO I=1,ND            !actualizamos XDD para futuras llamadas a FUNK's
            XDD(I)=XD(I)
          END DO
! -> REFINAMOS y --------------------------------------------------------------
        ELSEIF(CREF.EQ.'3')THEN
          WRITE(*,100)'Valor inicial  en Y: '
          WRITE(*,*) YD(NREF)
          XX0(1)=YD(NREF)                      !valores iniciales para DOWNHILL
          IF(YD(NREF).NE.0.0)THEN
            DXX0(1)=YD(NREF)*0.05
          ELSEIF(YD(1).NE.YD(ND))THEN
            DXX0(1)=(YD(1)-YD(ND))/5.
          ELSE
            DXX0(1)=1. !que remedio
          END IF
          WRITE(*,100)'Running DOWNHILL...'
          CALL DOWNHILL(1,XX0,DXX0,FUNKSPLFIT2,1.0,0.5,2.0,YRMSTOL,XX,DXX,NEVAL,500)
          WRITE(*,110)'      no. of function evaluations: ',NEVAL
          YD(NREF)=XX(1)
          SIGMA=SQRT(FUNKSPLFIT2(XX))
          WRITE(*,100)'Valor refinado en Y: '
          WRITE(*,*) YD(NREF)
! -> refinamos todos los nodos menos el primero y el ultimo --------------------
        ELSEIF(CREF.EQ.'W')THEN
          CALL RANSPL(ND,NRANND)            !ordenamos los Knots aleatoriamente
          NITER=NITER+1
          WRITE(*,111)'>>> ITERATION #',NITER
          WRITE(*,100)'     --> '
          DO I=1,ND-1                !mostramos el orden aleatorio de los Knots
            WRITE(CDUMMY,*)NRANND(I)
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,100)CDUMMY(1:L)//','
          END DO
          WRITE(CDUMMY,*)NRANND(ND)
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)CDUMMY(1:L)
          DO H=1,ND
            NREF=NRANND(H)
            IF(NREF.EQ.1)THEN
              YD(NREF)=YINI
            ELSEIF(NREF.EQ.ND)THEN
              YD(NREF)=YFIN
            ELSE                                            !refinamos en X e Y
              XX0(1)=XD(NREF)                  !valores iniciales para DOWNHILL
              IF(XD(NREF-1).NE.XD(NREF))THEN
                DXX0(1)=(XD(NREF)-XD(NREF-1))*0.05
              ELSEIF(XD(NREF+1).NE.XD(NREF))THEN
                DXX0(1)=(XD(NREF)-XD(NREF+1))*0.05
              ELSE
                DXX0(1)=1. !que remedio
              END IF
              XX0(2)=YD(NREF)
              IF(YD(NREF).NE.0.0)THEN
                DXX0(2)=YD(NREF)*0.05
              ELSEIF(YD(NREF-1).NE.YD(NREF))THEN
                DXX0(2)=(YD(NREF)-YD(NREF-1))*0.2
              ELSEIF(YD(NREF+1).NE.YD(NREF))THEN
                DXX0(2)=(YD(NREF)-YD(NREF+1))*0.2
              ELSE
                DXX0(2)=1. !que remedio
              END IF
              WRITE(*,100)'Running DOWNHILL...'
              CALL DOWNHILL(2,XX0,DXX0,FUNKSPLFIT3,1.0,0.5,2.0,YRMSTOL,XX,DXX,NEVAL,500)
              WRITE(*,110)'      no. of function evaluations: ',NEVAL
              XD(NREF)=XX(1)
              YD(NREF)=XX(2)
              SIGMA=SQRT(FUNKSPLFIT3(XX))
              DO I=1,ND        !actualizamos XDD para futuras llamadas a FUNK's
                XDD(I)=XD(I)
              END DO
            END IF
          END DO
        END IF
!------------------------------------------------------------------------------
        GOTO 20
!------------------------------------------------------------------------------
900     DO K=1,NOUT
          YOUT(K)=0.
        END DO
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
111     FORMAT(A,I6,$)
        END
!
!******************************************************************************
! Funcion para minimizar las coordenadas Y de todos los Knots
        REAL FUNCTION FUNKSPLFIT(X)
        IMPLICIT NONE
        INCLUDE 'largest.inc'
        INTEGER NDMAX
        PARAMETER (NDMAX=20)
        REAL X(NDMAX)
!
        INTEGER I
        INTEGER NF,ND
        INTEGER I0
        REAL XF(NXYMAX),YF(NXYMAX),YF0
        REAL XDD(NDMAX)
        REAL S(NDMAX),A(NDMAX),B(NDMAX),C(NDMAX)
        DOUBLE PRECISION F
!
        COMMON/BLKFUNK1/NF
        COMMON/BLKFUNK2/XF,YF
        COMMON/BLKFUNK3/ND
        COMMON/BLKFUNK4/XDD
!------------------------------------------------------------------------------
        S(1)=0.
        S(ND)=0.
        CALL CUBSPL(XDD,X,ND,4,S,A,B,C) !...............................IMODE=1
        F=0.D0
        I0=1                     !la primera vez busca en el inicio de la tabla
        DO I=1,NF
          CALL CUBSPLX(XDD,X,A,B,C,ND,I0,XF(I),YF0)
          F=F+(DBLE(YF0)-DBLE(YF(I)))*(DBLE(YF0)-DBLE(YF(I)))
        END DO
        F=F/DBLE(NF)
        FUNKSPLFIT=REAL(F)
!
        END
!
!******************************************************************************
! Funcion para minimizar la coordenada X de un solo Knot
        REAL FUNCTION FUNKSPLFIT1(X)
        IMPLICIT NONE
        INCLUDE 'largest.inc'
        INTEGER NDMAX
        PARAMETER (NDMAX=20)
        REAL X(NDMAX)
!
        INTEGER I
        INTEGER I0
        INTEGER NF,ND
        INTEGER NREF
        REAL XF(NXYMAX),YF(NXYMAX),YF0
        REAL XDD(NDMAX)
        REAL S(NDMAX),A(NDMAX),B(NDMAX),C(NDMAX)
        REAL YD(NDMAX)
        DOUBLE PRECISION F
!
        COMMON/BLKFUNK1/NF
        COMMON/BLKFUNK2/XF,YF
        COMMON/BLKFUNK3/ND
        COMMON/BLKFUNK4/XDD
        COMMON/BLKFUNK5/YD
        COMMON/BLKFUNK6/NREF
!------------------------------------------------------------------------------
        IF(X(1).LE.XDD(NREF-1)) GOTO 900
        IF(X(1).GE.XDD(NREF+1)) GOTO 900
        XDD(NREF)=X(1)
        S(1)=0.
        S(ND)=0.
        CALL CUBSPL(XDD,YD,ND,4,S,A,B,C) !..............................IMODE=1
        F=0.D0
        I0=1                     !la primera vez busca en el inicio de la tabla
        DO I=1,NF
          CALL CUBSPLX(XDD,YD,A,B,C,ND,I0,XF(I),YF0)
          F=F+(DBLE(YF0)-DBLE(YF(I)))*(DBLE(YF0)-DBLE(YF(I)))
        END DO
        F=F/DBLE(NF)
        FUNKSPLFIT1=REAL(F)
        RETURN
!
900     FUNKSPLFIT1=1.E10
        END
!
!******************************************************************************
! Funcion para minimizar la coordenada Y de un solo Knot
        REAL FUNCTION FUNKSPLFIT2(X)
        IMPLICIT NONE
        INCLUDE 'largest.inc'
        INTEGER NDMAX
        PARAMETER (NDMAX=20)
        REAL X(NDMAX)
!
        INTEGER I
        INTEGER I0
        INTEGER NF,ND
        INTEGER NREF
        REAL XF(NXYMAX),YF(NXYMAX),YF0
        REAL XDD(NDMAX)
        REAL S(NDMAX),A(NDMAX),B(NDMAX),C(NDMAX)
        REAL YD(NDMAX)
        DOUBLE PRECISION F
!
        COMMON/BLKFUNK1/NF
        COMMON/BLKFUNK2/XF,YF
        COMMON/BLKFUNK3/ND
        COMMON/BLKFUNK4/XDD
        COMMON/BLKFUNK5/YD
        COMMON/BLKFUNK6/NREF
!------------------------------------------------------------------------------
        YD(NREF)=X(1)
        S(1)=0.
        S(ND)=0.
        CALL CUBSPL(XDD,YD,ND,4,S,A,B,C) !..............................IMODE=1
        F=0.D0
        I0=1                     !la primera vez busca en el inicio de la tabla
        DO I=1,NF
          CALL CUBSPLX(XDD,YD,A,B,C,ND,I0,XF(I),YF0)
          F=F+(DBLE(YF0)-DBLE(YF(I)))*(DBLE(YF0)-DBLE(YF(I)))
        END DO
        F=F/DBLE(NF)
        FUNKSPLFIT2=REAL(F)
!
        END
!
!******************************************************************************
! Funcion para minimizar las coordenadas X e Y de un solo Knot
        REAL FUNCTION FUNKSPLFIT3(X)
        IMPLICIT NONE
        INCLUDE 'largest.inc'
        INTEGER NDMAX
        PARAMETER (NDMAX=20)
        REAL X(NDMAX)
!
        INTEGER I
        INTEGER I0
        INTEGER NF,ND
        INTEGER NREF
        REAL XF(NXYMAX),YF(NXYMAX),YF0
        REAL XDD(NDMAX)
        REAL S(NDMAX),A(NDMAX),B(NDMAX),C(NDMAX)
        REAL YD(NDMAX)
        DOUBLE PRECISION F
!
        COMMON/BLKFUNK1/NF
        COMMON/BLKFUNK2/XF,YF
        COMMON/BLKFUNK3/ND
        COMMON/BLKFUNK4/XDD
        COMMON/BLKFUNK5/YD
        COMMON/BLKFUNK6/NREF
!------------------------------------------------------------------------------
        IF(X(1).LE.XDD(NREF-1)) GOTO 900
        IF(X(1).GE.XDD(NREF+1)) GOTO 900
        XDD(NREF)=X(1)
        YD(NREF)=X(2)
        S(1)=0.
        S(ND)=0.
        CALL CUBSPL(XDD,YD,ND,4,S,A,B,C) !..............................IMODE=1
        F=0.D0
        I0=1                     !la primera vez busca en el inicio de la tabla
        DO I=1,NF
          CALL CUBSPLX(XDD,YD,A,B,C,ND,I0,XF(I),YF0)
          F=F+(DBLE(YF0)-DBLE(YF(I)))*(DBLE(YF0)-DBLE(YF(I)))
        END DO
        F=F/DBLE(NF)
        FUNKSPLFIT3=REAL(F)
        RETURN
!
900     FUNKSPLFIT3=1.E10
        END
!
!******************************************************************************
! Ordena aleatoriamente los numeros 1,2,...N en la matriz IX(). Esta subrutina
! es util para optimizar los Knots en orden aleatorio para evitar efectos
! sistematicos
        SUBROUTINE RANSPL(N,IX)
        IMPLICIT NONE
        INTEGER NMAX
        PARAMETER(NMAX=100)
        INTEGER N
        INTEGER IX(N),IS(NMAX)
!
        INTEGER NSEED
        INTEGER I,K,M,L
        REAL RANDOMNUMBER
        COMMON/BLKNSEED/NSEED
!------------------------------------------------------------------------------
        IF(N.GT.NMAX)THEN
          WRITE(*,101)'FATAL ERROR: N.GT.NMAX in RANSPL'
        END IF
!
        DO I=1,N
          IS(I)=I
        END DO
!
        K=N
        I=0
!
        DO WHILE(K.GE.1)
          M=INT(RANDOMNUMBER(NSEED)*REAL(K))+1    !numero aleatorio entre 1 y K
          I=I+1
          IX(I)=IS(M)
          IF(M.LT.K)THEN
            DO L=M,K-1
              IS(L)=IS(L+1)
            END DO
          END IF
          K=K-1
        END DO
!
101     FORMAT(A)
        END
