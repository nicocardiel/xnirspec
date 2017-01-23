C------------------------------------------------------------------------------
C Version 10-September-1999                                      File: polfit.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE  POLFIT(X,Y,SIGMAY,NPTS,NTERMS,MODE,A,CHISQR)
C
C >>> This subroutine is based on the subroutine from Data Reduction and 
C Error Analysis for the Physical Sciences (Bevington, 1969) <<<
C
C        LEAST-SQUARES FIT TO A POLYNOMIAL
C        INPUT: X  -  ARRAY FOR INDEPENDENT VARIABLE
C               Y  -  ARRAY FOR DEPENDENT VARIABLE
C               SIGMAY  -  STANDARD DEVIATIONS FOR Y DATA POINTS
C               NPTS  -  NUMBER OF PAIRS OF DATA POINTS
C               NTERMS  - NUMBER OF COEFFICIENTS (DEGREE + 1)
C               MODE  -  METHOD OF WEIGHTING (0 = NO WEIGHTING)
C        OUTPUT:A  - ARRAY OF COEFFICIENTS
C               CHISQR  -  REDUCED CHI SQUARE FOR FIT
C
C        IT USES FUNCTION DETERM TO EVALUATE DETERMINANT OF MATRIX
C        SUPPORTS NTERM UP TO 20
C        FOR DETAILS SEE BEVINGTON(1969)
C
Comment
C------------------------------------------------------------------------------
C Nota: en la llamada a esta subrutina no introducir un "0" en CHISQR aunque
C no nos interese el valor de esta variable a la salida de la subrutina. En
C caso contrario obtenemos resultamos impredecibles en los ajustes.
        SUBROUTINE POLFIT(X,Y,SIGMAY,NPTS,NTERMS,MODE,A,CHISQR)
        implicit none
        integer npts,nterms,mode
        real chisqr
C
        integer i,n,nmax,j,k,l
        real xi,yi,weight,free
        DOUBLE PRECISION SUMX(39),SUMY(20)
        DOUBLE PRECISION XTERM,YTERM,CHISQ
        real cx1,cx2,cy1,cy2
        real xmin,xmax,ymin,ymax
        DOUBLE PRECISION DETERM,DELTA,ARRAY(20,20)
        REAL X(NPTS),Y(NPTS),SIGMAY(NPTS),A(NTERMS)
        real aa(20)
        double precision combpf
C------------------------------------------------------------------------------
        cx1=1. !evita WARNING de compilacion
        cx2=0. !evita WARNING de compilacion
        cy1=1. !evita WARNING de compilacion
        cy2=0. !evita WARNING de compilacion
c normalizamos X,Y para que tengan valores en el intervalo [-1,+1], solo
c cuando MODE=0
        if(mode.eq.0)then
          xmax=x(1)
          ymax=y(1)
          xmin=x(1)
          ymin=y(1)
          if(npts.gt.1)then
            do i=2,npts
              if(x(i).gt.xmax) xmax=x(i)
              if(y(i).gt.ymax) ymax=y(i)
              if(x(i).lt.xmin) xmin=x(i)
              if(y(i).lt.ymin) ymin=y(i)
            end do
          end if
          if(xmin.eq.xmax)then
            cx1=1.
            cx2=0.
          else
            cx1=2./(xmax-xmin)
            cx2=(xmax+xmin)/(xmax-xmin)
          end if
          !como medida de precaucion, evitamos el caso xmax=xmin
          if(cx1.eq.0.0)then
            cx1=1.0
            cx2=0.0
          end if
          do i=1,npts
            x(i)=x(i)*cx1-cx2
          end do
          if(ymin.eq.ymax)then
            cy1=1.
            cy2=0.
          else
            cy1=2./(ymax-ymin)
            cy2=(ymax+ymin)/(ymax-ymin)
          end if
          !como medida de precaucion, evitamos el caso ymax=ymin
          if(cy1.eq.0.0)then
            cy1=1.0
            cy2=0.0
          end if
          do i=1,npts
            y(i)=y(i)*cy1-cy2
          end do
        end if
C------------------------------------------------------------------------------
C        ACCUMULATE WEIGHTED SUMS
C
        NMAX=2*NTERMS-1
        DO N=1,NMAX
          SUMX(N)=0.d0
        END DO
        DO J=1,NTERMS
          SUMY(J)=0.d0
        END DO
        CHISQ=0.
        DO I=1,NPTS
          XI=X(I)
          YI=Y(I)
          IF (MODE) 32,37,39
32        IF (YI) 35,37,33
33        WEIGHT=1./YI
          GO TO 41
35        WEIGHT=1./(-YI)
          GO TO 41
37        WEIGHT=1.
          GO TO 41
39        WEIGHT=1./SIGMAY(I)**2
41        XTERM=dble(WEIGHT)
          DO N=1,NMAX
            SUMX(N)=SUMX(N)+XTERM
            XTERM=XTERM*dble(XI)
          END DO
          YTERM=dble(WEIGHT)*dble(YI)
          DO N=1,NTERMS
            SUMY(N)=SUMY(N)+YTERM
            YTERM=YTERM*dble(XI)
          END DO
          CHISQ=CHISQ+dble(WEIGHT)*dble(YI)**2.d0
        END DO
C
C        CONSTRUCT MATRICES AND CALCULATE COEFFICIENTS
C
        DO J=1,NTERMS
          DO K=1,NTERMS
            N=J+K-1
            ARRAY(J,K)=SUMX(N)
          END DO
        END DO
        DELTA=DETERM(ARRAY,NTERMS)
        IF (DELTA) 61,57,61
57      CHISQR=0.
        DO J=1,NTERMS
          A(J)=0.
        END DO
        GO TO 80
61      DO L=1,NTERMS
          DO J=1,NTERMS
            DO K=1,NTERMS
              N=J+K-1
              ARRAY(J,K)=SUMX(N)
            END DO
            ARRAY(J,L)=SUMY(J)
          END DO
          A(L)=real(DETERM(ARRAY,NTERMS)/DELTA)
        END DO
C
C        CALCULATE CHI SQUARE
C
        DO J=1,NTERMS
          CHISQ=CHISQ-2.d0*dble(A(J))*SUMY(J)
          DO K=1,NTERMS
            N=J+K-1
            CHISQ=CHISQ+dble(A(J))*dble(A(K))*SUMX(N)
          END DO
        END DO
C CHISQR, modificado por NCL cuando NPTS=NTERMS
        if(npts-nterms.gt.0)then
          free=real(npts-nterms)
          CHISQR=real(CHISQ/dble(FREE))
        else
          chisqr=0.
        end if
80      continue
c------------------------------------------------------------------------------
c deshacemos el cambio de variable
        if(mode.eq.0)then
          chisqr=0.       !porque el cambio de variable tambien afecta a CHISQR
          do i=1,npts
            x(i)=(x(i)+cx2)/cx1
            y(i)=(y(i)+cy2)/cy1
          end do
          if(cx2.eq.0.0)then
            do k=0,nterms-1
              a(k+1)=a(k+1)*(cx1**(k))
            end do
          else
            do k=0,nterms-1
              aa(k+1)=a(k+1)
            end do
            do k=0,nterms-1
              a(k+1)=0.
              do i=k,nterms-1
                a(k+1)=a(k+1)+aa(i+1)*real(combpf(i,i-k))*(cx1**k)*
     +           ((-cx2)**(i-k))
              end do
            end do
          end if
            a(1)=a(1)+cy2
           do k=0,nterms-1
            a(k+1)=a(k+1)/cy1
          end do
        end if
c
        END
C
C******************************************************************************
C******************************************************************************
C
        DOUBLE PRECISION FUNCTION DETERM(ARRAY,NORDER)
C
        DOUBLE PRECISION ARRAY(20,20),SAVE
C
        DETERM=1.d0
        DO 50 K=1,NORDER
          IF(ARRAY(K,K)) 41,21,41
21        DO 23 J=K,NORDER
ccccc            IF(ARRAY(J,K)) 31,23,31  !HORROR!!!!!!!!!!!!!!!!!!!!!!!!
            IF(ARRAY(K,J)) 31,23,31
23        CONTINUE
          DETERM=0.
          GO TO 60
31        DO I=K,NORDER
            SAVE=ARRAY(I,J)
            ARRAY(I,J)=ARRAY(I,K)
            ARRAY(I,K)=SAVE
          END DO
          DETERM=-DETERM
41        DETERM=DETERM*ARRAY(K,K)
          IF (K-NORDER) 43,50,50
43        K1=K+1
          DO I=K1,NORDER
            DO J=K1,NORDER
              ARRAY(I,J)=ARRAY(I,J)-ARRAY(I,K)*ARRAY(K,J)/ARRAY(K,K)
            END DO
          END DO
50      CONTINUE
60      RETURN
        END
