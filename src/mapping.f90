        SUBROUTINE MAPPING(LBOUNDARY,LMAPPING)
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_COEFFBL
        USE Dynamic_Array_COEFFBA
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_coeffbl.inc'
        INCLUDE 'interface_coeffba.inc'
! subroutine arguments
        LOGICAL LBOUNDARY,LMAPPING
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER READI,READILIM
        INTEGER TRUEBEG,TRUELEN
        INTEGER SYSTEMFUNCTION
        REAL READF
        REAL ARCLENGTH
        REAL FMEAN0,FMEAN2
        CHARACTER*255 READC
!
        INTEGER NMAPMAX
        PARAMETER (NMAPMAX=3) !second-degree bivariate approximation
        INTEGER NECUAMAX
        PARAMETER (NECUAMAX=(NMAPMAX+1)*NMAPMAX)
!
!
        INTEGER ISYSTEM
        INTEGER L1,L2
        INTEGER N,NECUA
        INTEGER I,J,L,LL,M,II,JJ,K,I1,I2,J1,J2
        INTEGER NFIT,NMAP
        INTEGER NP
        INTEGER NDEGBL(NXYMAX),NDEGBA(NXYMAX)  !polynomial degrees for boundary
        INTEGER NDEGBL00,NDEGBA00
        INTEGER NLINBL,NLINBA                  !number of arc lines in boundary
        INTEGER NREFBL,NREFBA
        INTEGER ORDER(NECUAMAX),IOK,IPAR
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        INTEGER NDELTAX,NDELTAY
        INTEGER NWIDTH
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NCBUFF
        INTEGER NEXTINFO
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
!delete REAL COEFFBL(20,NXYMAX),COEFFBA(20,NXYMAX)      !pol. coef. in boundary
        REAL X(NXYMAX),Y(NXYMAX),U(NXYMAX),V(NXYMAX)
        REAL XX,YY,UU,VV
        REAL X0,Y0,X1,Y1
        REAL XX0(1),YY0(1)
        REAL SX,SY
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XP(NXYMAX),YP(NXYMAX),RP(NXYMAX)
        REAL XP_(NXYMAX),YP_(NXYMAX)
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL XMINBL(NXYMAX),XMAXBL(NXYMAX)
        REAL YMINBL(NXYMAX),YMAXBL(NXYMAX)
        REAL XMINBA(NXYMAX),XMAXBA(NXYMAX)
        REAL YMINBA(NXYMAX),YMAXBA(NXYMAX)
        REAL FSCALE
        REAL SUM
        REAL AIJ(NECUAMAX),BIJ(NECUAMAX)
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
        REAL A(NECUAMAX,NECUAMAX),B(NECUAMAX)
        REAL SCALEROW(NECUAMAX)
        REAL CHISQR1,CHISQR2
        REAL XC,YC
        REAL ALPHA0,ALPHA1,ALPHA2
        REAL ALPHA0_,ALPHA1_,ALPHA2_
        REAL BETA(5) !coeficientes para polinomio de cuarto grado
        REAL XMIN_,XMAX_,YMIN_,YMAX_
        REAL AMP,SIGMA,FACTOR
        REAL EEX0,EEY0,EEAMP,EESIGMA
        REAL FMEAN,FSIGMA
        CHARACTER*1 CEXT,COPC,CH,CMORE
        CHARACTER*50 CDUMMY
        CHARACTER*255 FILENAME
        LOGICAL LBOUNDL,LBOUNDA
        LOGICAL LOGFILE,LOGFILERR
!
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKIMAGEN2/NCBUFF
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKBOUND1/NLINBL,NDEGBL,NDEGBL00
        COMMON/BLKBOUND1B/NLINBA,NDEGBA,NDEGBA00
!delete COMMON/BLKBOUND2/COEFFBL
!delete COMMON/BLKBOUND2B/COEFFBA
        COMMON/BLKBOUND3/LBOUNDL,LBOUNDA
        COMMON/BLKBOUND4/XMINBL,YMINBL,XMAXBL,YMAXBL
        COMMON/BLKBOUND5/XMINBA,YMINBA,XMAXBA,YMAXBA
        COMMON/BLKMAPPING0/NMAP
        COMMON/BLKMAPPING1/AIJ,BIJ
        COMMON/BLKMAPPING2/AIJ_,BIJ_
        COMMON/BLKMAPPING3/SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        COMMON/BLKMAPPING4/SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        COMMON/BLKFSCALE/FSCALE
        COMMON/BLKDEFAULTS2/NWIDTH
        COMMON/BLKPLIMITS/XMIN_,XMAX_,YMIN_,YMAX_
!------------------------------------------------------------------------------
        WRITE(*,101)
        WRITE(*,101) '(l) load mapping from file'
        IF(LMAPPING)THEN
          WRITE(*,101) '(s) save mapping into file'
          WRITE(*,101) '(t) test current mapping'
        END IF
        IF(LBOUNDARY)THEN
          WRITE(*,101) '(c) compute mapping from current boundary'
        END IF
        WRITE(*,101) '(x) exit'
        IF(LMAPPING.AND.LBOUNDARY)THEN
          COPC(1:1)=READC('Option (l/s/t/c/x)','x','lstcx')
        ELSEIF(LMAPPING)THEN
          COPC(1:1)=READC('Option (l/s/t/x)','x','lstx')
        ELSEIF(LBOUNDARY)THEN
          COPC(1:1)=READC('Option (l/c/x)','x','lcx')
        ELSE
          COPC(1:1)=READC('Option (l/x)','x','lx')
        END IF
!------------------------------------------------------------------------------
        IF(COPC.EQ.'x')THEN
          RETURN
!------------------------------------------------------------------------------
        ELSEIF(COPC.EQ.'l')THEN
          LOGFILE=.FALSE.
          LOGFILERR=.FALSE.
          DO WHILE(.NOT.LOGFILE)
            FILENAME=READC('Mapping file name (none=EXIT)','*.mapp','@')
            IF((INDEX(FILENAME,'*').NE.0).OR.(INDEX(FILENAME,'?').NE.0))THEN
              L1=TRUEBEG(FILENAME)
              L2=TRUELEN(FILENAME)
              ISYSTEM=SYSTEMFUNCTION('ls '//FILENAME(L1:L2))
            ELSEIF(FILENAME.EQ.'none')THEN
              LOGFILE=.TRUE.
              LOGFILERR=.TRUE.
            ELSE
              INQUIRE(FILE=FILENAME,EXIST=LOGFILE)
              IF(.NOT.LOGFILE)THEN
                WRITE(*,101) '***ERROR***'
                WRITE(*,100) '=> This file does not exist.'
                WRITE(*,100) ' Try again. (press <CR>...)'
                READ(*,*)
              END IF
            END IF
          END DO
          IF(.NOT.LOGFILERR)THEN                   !the file name is not 'none'
            WRITE(*,100) 'Reading file...'
            OPEN(10,FILE=FILENAME,STATUS='OLD',FORM='FORMATTED')
            READ(10,*) !skip line
            READ(10,*) FSCALE
            FSCALE=1./FSCALE
            READ(10,*) !skip line
            READ(10,*) NMAP
            NECUA=(NMAP+1)*NMAP
            READ(10,*) !skip line
            DO I=1,NECUA
              READ(10,*) AIJ(I)
            END DO
            READ(10,*) !skip line
            DO I=1,NECUA
              READ(10,*) BIJ(I)
            END DO
            READ(10,*) !skip line
            DO I=1,NECUA
              READ(10,*) AIJ_(I)
            END DO
            READ(10,*) !skip line
            DO I=1,NECUA
              READ(10,*) BIJ_(I)
            END DO
            READ(10,*) !skip line
            READ(10,*) SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
            READ(10,*) !skip line
            READ(10,*) SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
            CLOSE(10)
            WRITE(*,101) ' ...OK!'
            I1=SYMINGRID-SYMINEXTG !borde inferior del primer pixel
            I2=SYMAXGRID+SYMAXEXTG !borde superior del ultimo pixel
            J1=SXMINGRID-SXMINEXTG !borde izquierdo del primer pixel
            J2=SXMAXGRID+SXMAXEXTG !borde derecho del ultimo pixel
! OJO: el numero de pixels no es I2-I1+1, sino que es I2-I1; por eso los bucles
! van desde I1 hasta I2-1. I1,I2,J1,J2 se refieren a los bordes de los pixels,
! no al numero de pixel. Hay un offset de 0.5 entre el mapeado realizado con
! los coeficientes AIJ,BIJ,AIJ_,BIJ_ y el numero de pixel en la imagen
! original.
            NDELTAX=J2-J1
            NDELTAY=I2-I1
            WRITE(*,100) '> Image size X*Y: '
            WRITE(CDUMMY,'(I10,A1,I10)') NDELTAX,'*',NDELTAY
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            WRITE(*,101) CDUMMY(1:L)
            LMAPPING=.TRUE.
            CALL PMAPP
          END IF
          RETURN
!------------------------------------------------------------------------------
        ELSEIF(COPC.EQ.'s')THEN
          LOGFILE=.TRUE.
          LOGFILERR=.FALSE.
          DO WHILE(LOGFILE)
            FILENAME=READC('Mapping file name (none=EXIT)','*.mapp','@')
            IF((INDEX(FILENAME,'*').NE.0).OR.(INDEX(FILENAME,'?').NE.0))THEN
              L1=TRUEBEG(FILENAME)
              L2=TRUELEN(FILENAME)
              ISYSTEM=SYSTEMFUNCTION('ls '//FILENAME(L1:L2))
            ELSEIF(FILENAME.EQ.'none')THEN
              LOGFILE=.FALSE.
              LOGFILERR=.TRUE.
            ELSE
              INQUIRE(FILE=FILENAME,EXIST=LOGFILE)
              IF(LOGFILE)THEN
                WRITE(*,101) '***ERROR***'
                WRITE(*,100) '=> This file does already exist.'
                WRITE(*,100) ' Try again. (press <CR>...)'
                READ(*,*)
              END IF
            END IF
          END DO
          IF(.NOT.LOGFILERR)THEN                   !the file name is not 'none'
            WRITE(*,100) 'Writting file...'
            OPEN(10,FILE=FILENAME,STATUS='NEW',FORM='FORMATTED')
            WRITE(10,101) '#Scale factor FSCALE'
            WRITE(10,*) 1./FSCALE
            WRITE(10,101) '#Bivariate polynomial degree'
            WRITE(10,*) NMAP
            NECUA=(NMAP+1)*NMAP
            WRITE(10,101) '# aij:'
            DO I=1,NECUA
              WRITE(10,*) AIJ(I)
            END DO
            WRITE(10,101) '# bij:'
            DO I=1,NECUA
              WRITE(10,*) BIJ(I)
            END DO
            WRITE(10,101) '# aij_:'
            DO I=1,NECUA
              WRITE(10,*) AIJ_(I)
            END DO
            WRITE(10,101) '# bij_:'
            DO I=1,NECUA
              WRITE(10,*) BIJ_(I)
            END DO
            WRITE(10,101) '#Grid limits (no extrapolation)'
            WRITE(10,*) SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
            WRITE(10,101) '#Grid extrapolations'
            WRITE(10,*) SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
            CLOSE(10)
            WRITE(*,101) ' ...OK!'
          END IF
          RETURN
!------------------------------------------------------------------------------
        ELSEIF(COPC.EQ.'t')THEN
          IF(NMAP.NE.2)THEN
            WRITE(*,100) 'WARNING: this test only available when '
            WRITE(*,101) 'NMAP=2'
            WRITE(*,100) '(Press <CR> to continue...)'
            READ(*,*)
            RETURN
          END IF
          WRITE(CDUMMY,*) NWIDTH
          NWIDTH=READI('NWIDTH',CDUMMY)
          CMORE(1:1)=READC('Plot individual fits (y/n)','n','yn')
          CH='A'
          DO WHILE(CH.NE.'X')
            WRITE(*,101) 'Click mouse in sky line...'
            CALL PGSCI(5)
            CALL RPGBAND(7,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IF(CH.EQ.'X')THEN
              WRITE(*,101) '...operation canceled!'
            ELSE
              WRITE(*,101) '...OK!'
              CALL FMAP(NMAP,AIJ,BIJ,XC,YC,UU,VV)
              DO I=1,NXYMAX
                VV=REAL(SYMINGRID-SYMINEXTG)+REAL(I-1)/REAL(NXYMAX-1)*REAL(SYMAXGRID+SYMAXEXTG-SYMINGRID+SYMINEXTG)
                CALL FMAP(NMAP,AIJ_,BIJ_,UU,VV,XP(I),YP(I))
              END DO
              CALL PGSCI(6)
              CALL PGLINE(NXYMAX,XP,YP) !dibujamos posicion estimada
              CALL PGSCI(1)
              VV=REAL(SYMAXGRID+SYMAXEXTG) !extremo superior de la linea
              CALL FMAP(NMAP,AIJ_,BIJ_,UU,VV,XX,YY)
              I2=NINT(YY) !extremo superior de la linea
              VV=REAL(SYMINGRID-SYMINEXTG) !extremo inferior de la linea
              CALL FMAP(NMAP,AIJ_,BIJ_,UU,VV,XX,YY)
              I1=NINT(YY) !extremo inferior de la linea
              ALPHA0_=AIJ_(1)+AIJ_(4)*(UU*FSCALE)+AIJ_(6)*(UU*FSCALE)*(UU*FSCALE)
              ALPHA1_=AIJ_(2)+AIJ_(5)*(UU*FSCALE)
              ALPHA2_=AIJ_(3)
              NEXTINFO=0
              DO I=I1,I2
                YY=REAL(I)
                ALPHA0=BIJ(1)+BIJ(2)*(YY*FSCALE)+BIJ(3)*(YY*FSCALE)*(YY*FSCALE)
                ALPHA1=BIJ(4)+BIJ(5)*(YY*FSCALE)
                ALPHA2=BIJ(6)
                BETA(1)=ALPHA0_+ALPHA1_*ALPHA0+ALPHA2_*ALPHA0*ALPHA0
                BETA(2)=ALPHA1_*ALPHA1+2.*ALPHA2_*ALPHA0*ALPHA1-1.
                BETA(3)=ALPHA1_*ALPHA2+2.*ALPHA2_*ALPHA0*ALPHA2+ALPHA2_*ALPHA1*ALPHA1
                BETA(4)=2.*ALPHA2_*ALPHA1*ALPHA2
                BETA(5)=ALPHA2_*ALPHA2*ALPHA2
                X0=XX*FSCALE !valor inicial para Newton-Raphson
                CALL FPOLYINV(4,BETA,0.,X0,XX)
                XX=XX/FSCALE
                J1=NINT(XX)-NWIDTH/2
                J2=NINT(XX)+NWIDTH/2
                IF(J1.LT.1)THEN
                  J1=1
                  J2=J1+NWIDTH-1
                ELSEIF(J2.GT.NAXIS(1,NCBUFF))THEN
                  J2=NAXIS(1,NCBUFF)
                  J1=J2-NWIDTH+1
                END IF
                DO J=J1,J2
                  XF(J-J1+1)=REAL(J)
                  YF(J-J1+1)=IMAGEN(J,I,NCBUFF)
                END DO
                CALL GAUSCFIT(NWIDTH,XF,YF,X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,1.E-6)
                XP(I-I1+1)=X0
                YP(I-I1+1)=REAL(I)
                RP(I-I1+1)=XP(I-I1+1)-XX
                IF(CMORE.EQ.'y')THEN
                  XMIN_=REAL(J1)-0.6
                  XMAX_=REAL(J2)+0.6
                  CALL SUBPLOT(NWIDTH,1,NWIDTH,XF,YF,XF,YF,.FALSE.,.TRUE.,.FALSE.,.FALSE.,'y-axis','signal',' ',3,201,1.0)
                  DO II=1,NXYMAX
                    XP_(II)=REAL(J1)+REAL(II-1)/REAL(NXYMAX-1)*REAL(J2-J1)
                    FACTOR=(XP_(II)-X0)*(XP_(II)-X0)/(2.*SIGMA*SIGMA)
                    IF(FACTOR.GT.60.)THEN
                      YP_(II)=Y0
                    ELSE
                      YP_(II)=Y0+AMP*EXP(-FACTOR)
                    END IF
                  END DO
                  CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,XP_,YP_,XP_,YP_,.FALSE.,.FALSE.,5,101,1.0)
                  print*,'scans: ',i,'(',i1,',',i2,')'
                  print*,'x0...: ',x0,eex0
                  print*,'sigma: ',sigma,eesigma
                  print*,'amp..: ',amp,eeamp
                  print*,'y0...: ',y0,eey0
                  IF(CMORE.EQ.'y')THEN
                    CMORE(1:1)=READC('More plots (y/n)','y','yn')
                  END IF
                END IF
                CALL SHOWPERC(I1,I2,1,I,NEXTINFO)
              END DO
              NP=I2-I1+1
              CALL PGSCI(4)
              CALL PGPOINT(NP,XP,YP,1)
              CALL PGSCI(1)
              FMEAN=FMEAN0(NP,RP,FSIGMA)
              WRITE(*,100) '> Mean, r.m.s.............: '
              WRITE(*,*) FMEAN,FSIGMA
              FMEAN=FMEAN2(NP,RP,3.0,FSIGMA)
              WRITE(*,100) '> Mean, r.m.s. (< 3 sigma): '
              WRITE(*,*) FMEAN,FSIGMA
              YMIN_=FMEAN-5.*FSIGMA
              YMAX_=FMEAN+5.*FSIGMA
              CALL SUBPLOT(NP,1,NP,YP,RP,YP,RP,.TRUE.,.FALSE.,.FALSE.,.FALSE.,'spatial direction','offset residuals',' ',3,1,1.5)
            END IF
          END DO
          RETURN
        END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Aqui empieza la opcion 'c': calcular el mapping a partir del boundary
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! dibujamos el boundary con todas las intersecciones
        CALL PBOUND
!..............................................................................
! origen de medida de longitudes de arco
        NREFBL=READILIM('NREFBL','1',1,NLINBL)
        NREFBA=READILIM('NREFBA','1',1,NLINBA)
        CALL INTERSEC(NDEGBA(NREFBA),COEFFBA(1,NREFBA),NDEGBL(NREFBL),COEFFBL(1,NREFBL),0.5*(XMINBA(NREFBA)+XMAXBA(NREFBA)),X0,Y0)
        CALL PGSCI(3)
        !usamos un array unidimensional porque el compilador
        !gfortran-mp-10 da error al usar un escalar en lugar
        !de una matriz
        XX0(1)=X0
        YY0(1)=Y0
        CALL PGPOINT(1,XX0,YY0,24)
!..............................................................................
! calculamos todas las intersecciones y almacenamos los puntos para el ajuste
        CALL PGSCI(7)
        K=0
        DO L=1,NLINBA
          DO LL=1,NLINBL
            K=K+1
!
            IF(K.GT.NXYMAX)THEN
              WRITE(*,101) '***FATAL ERROR***'
              WRITE(*,101) '=> No. of points for fit > NXYMAX'
              WRITE(*,100) '=> NLINBL,NLINBA: '
              WRITE(*,*) NLINBL,NLINBA
              WRITE(*,101) '=> Redim X0 and Y0'
              WRITE(*,100) '(press <CR> to continue...)'
              READ(*,*)
              RETURN
            END IF
!
            CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),NDEGBL(LL),COEFFBL(1,LL),0.5*(XMINBA(L)+XMAXBA(L)),X(K),Y(K))
            CALL PGPOINT(1,X(K),Y(K),17)
!
            CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),NDEGBL(NREFBL),COEFFBL(1,NREFBL),0.5*(XMINBA(L)+XMAXBA(L)),X1,Y1)
            U(K)=ARCLENGTH(NDEGBL(NREFBL),COEFFBL(1,NREFBL),X0,X1,100)
!
            CALL INTERSEC(NDEGBA(NREFBA),COEFFBA(1,NREFBA),NDEGBL(LL),COEFFBL(1,LL),0.5*(XMINBA(NREFBA)+XMAXBA(NREFBA)),X1,Y1)
            V(K)=ARCLENGTH(NDEGBA(NREFBA),COEFFBA(1,NREFBA),Y0,Y1,100)
          END DO
        END DO
        NFIT=K
        CALL PGSCI(1)
!..............................................................................
! normalizamos las variables
        WRITE(CDUMMY,*) NXYMAX
        FSCALE=READF('FSCALE',CDUMMY)
        FSCALE=1./FSCALE
        DO K=1,NFIT
          X(K)=X(K)*FSCALE
          Y(K)=Y(K)*FSCALE
          U(K)=U(K)*FSCALE
          V(K)=V(K)*FSCALE
        END DO
!..............................................................................
        N=READILIM('Bivariate approximation degree','2',0,NMAPMAX)
        NECUA=(N+1)*N
        NMAP=N
!..............................................................................
! matriz del primer sistema de ecuaciones a resolver
! II: numero de ecuacion
! JJ: numero de incognita dentro de la ecuacion
! A(II,JJ): coeficiente de la ecuacion II e incognita JJ
        II=0
        DO L=0,N
          DO M=0,N-L
            II=II+1
            JJ=0
            DO I=0,N
              DO J=0,N-I
                JJ=JJ+1
                SUM=0.
                DO K=1,NFIT
                  SUM=SUM+(X(K)**(I+L))*(Y(K)**(J+M))
                END DO
                A(II,JJ)=SUM
              END DO
            END DO
            SUM=0.
            DO K=1,NFIT
              SUM=SUM+(X(K)**L)*(Y(K)**M)*U(K)
            END DO
            B(II)=SUM
          END DO
        END DO
! resolvemos el sistema de ecuaciones
        WRITE(*,100)'---> LU Descomposition...'
        CALL LUDCMP(A,NECUA,NECUAMAX,ORDER,SCALEROW,IOK,IPAR)
        WRITE(*,101)' OK!'
        WRITE(*,100)'---> Forward substitution and back substitution...'
        CALL LUSOLV(A,NECUA,NECUAMAX,ORDER,SCALEROW,B,AIJ)
        WRITE(*,101)' OK!'
!..............................................................................
! segundo sistema de ecuaciones a resolver (notar que la matriz
! del sistema de ecuaciones es la misma)
        II=0
        DO L=0,N
          DO M=0,N-L
            II=II+1
            SUM=0.
            DO K=1,NFIT
              SUM=SUM+(X(K)**L)*(Y(K)**M)*V(K)
            END DO
            B(II)=SUM
          END DO
        END DO
! resolvemos el sistema de ecuaciones
        WRITE(*,100)'---> Forward substitution and back substitution...'
        CALL LUSOLV(A,NECUA,NECUAMAX,ORDER,SCALEROW,B,BIJ)
        WRITE(*,101)' OK!'
!..............................................................................
! chisqr
        CHISQR1=0.
        CHISQR2=0.
        DO K=1,NFIT
          CALL FMAP(N,AIJ,BIJ,X(K)/FSCALE,Y(K)/FSCALE,UU,VV)
          CHISQR1=CHISQR1+(UU*FSCALE-U(K))*(UU*FSCALE-U(K))
          CHISQR2=CHISQR2+(VV*FSCALE-V(K))*(VV*FSCALE-V(K))
        END DO
!..............................................................................
! mostramos coeficientes aij
        II=0
        DO L=0,N
          DO M=0,N-L
            II=II+1
            WRITE(*,'(A10,I1,A1,I1,A3,$)') 'direct  a(',L,',',M,'): '
            WRITE(*,*) AIJ(II)
          END DO
        END DO
        WRITE(*,100) 'Residual standard deviation: '
        IF(NFIT.GT.NECUA)THEN
          WRITE(*,*) SQRT(CHISQR1/REAL(NFIT-NECUA))
        ELSE
          WRITE(*,*) 0.0
        END IF
! mostramos coeficientes bij
        II=0
        DO L=0,N
          DO M=0,N-L
            II=II+1
            WRITE(*,'(A10,I1,A1,I1,A3,$)') 'direct  b(',L,',',M,'): '
            WRITE(*,*) BIJ(II)
          END DO
        END DO
        WRITE(*,100) 'Residual standard deviation: '
        IF(NFIT.GT.NECUA)THEN
          WRITE(*,*) SQRT(CHISQR2/REAL(NFIT-NECUA))
        ELSE
          WRITE(*,*) 0.0
        END IF
!..............................................................................
! matriz del tercer sistema de ecuaciones a resolver
        II=0
        DO L=0,N
          DO M=0,N-L
            II=II+1
            JJ=0
            DO I=0,N
              DO J=0,N-I
                JJ=JJ+1
                SUM=0.
                DO K=1,NFIT
                  SUM=SUM+(U(K)**(I+L))*(V(K)**(J+M))
                END DO
                A(II,JJ)=SUM
              END DO
            END DO
            SUM=0.
            DO K=1,NFIT
              SUM=SUM+(U(K)**L)*(V(K)**M)*X(K)
            END DO
            B(II)=SUM
          END DO
        END DO
! resolvemos el sistema de ecuaciones
        WRITE(*,100)'---> LU Descomposition...'
        CALL LUDCMP(A,NECUA,NECUAMAX,ORDER,SCALEROW,IOK,IPAR)
        WRITE(*,101)' OK!'
        WRITE(*,100)'---> Forward substitution and back substitution...'
        CALL LUSOLV(A,NECUA,NECUAMAX,ORDER,SCALEROW,B,AIJ_)
        WRITE(*,101)' OK!'
!..............................................................................
! cuarto sistema de ecuaciones a resolver (notar que la matriz
! del sistema de ecuaciones es la misma)
        II=0
        DO L=0,N
          DO M=0,N-L
            II=II+1
            SUM=0.
            DO K=1,NFIT
              SUM=SUM+(U(K)**L)*(V(K)**M)*Y(K)
            END DO
            B(II)=SUM
          END DO
        END DO
! resolvemos el sistema de ecuaciones
        WRITE(*,100)'---> Forward substitution and back substitution...'
        CALL LUSOLV(A,NECUA,NECUAMAX,ORDER,SCALEROW,B,BIJ_)
        WRITE(*,101)' OK!'
!..............................................................................
! chisqr
        CHISQR1=0.
        CHISQR2=0.
        DO K=1,NFIT
          CALL FMAP(N,AIJ_,BIJ_,U(K)/FSCALE,V(K)/FSCALE,XX,YY)
          CHISQR1=CHISQR1+(XX*FSCALE-X(K))*(XX*FSCALE-X(K))
          CHISQR2=CHISQR2+(YY*FSCALE-Y(K))*(YY*FSCALE-Y(K))
        END DO
!..............................................................................
! mostramos coeficientes aij
        II=0
        DO L=0,N
          DO M=0,N-L
            II=II+1
            WRITE(*,'(A10,I1,A1,I1,A3,$)') 'inverse a(',L,',',M,'): '
            WRITE(*,*) AIJ_(II)
          END DO
        END DO
        WRITE(*,100) 'Residual standard deviation: '
        IF(NFIT.GT.NECUA)THEN
          WRITE(*,*) SQRT(CHISQR1/REAL(NFIT-NECUA))
        ELSE
          WRITE(*,*) 0.0
        END IF
! mostramos coeficientes bij
        II=0
        DO L=0,N
          DO M=0,N-L
            II=II+1
            WRITE(*,'(A10,I1,A1,I1,A3,$)') 'inverse b(',L,',',M,'): '
            WRITE(*,*) BIJ_(II)
          END DO
        END DO
        WRITE(*,100) 'Residual standard deviation: '
        IF(NFIT.GT.NECUA)THEN
          WRITE(*,*) SQRT(CHISQR2/REAL(NFIT-NECUA))
        ELSE
          WRITE(*,*) 0.0
        END IF
!..............................................................................
! dibujamos el grid con el nuevo ajuste
        CALL PGBBUF
        XMIN=U(NREFBL)/FSCALE
        XMAX=U((NLINBA-1)*NLINBL+NREFBL)/FSCALE
        DO LL=1,NLINBL
          DO I=1,NXYMAX
            SX=XMIN+REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)
            SY=V((NREFBA-1)*NLINBL+LL)/FSCALE
            CALL FMAP(N,AIJ_,BIJ_,SX,SY,XP(I),YP(I))
          END DO
          CALL PGSCI(5)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END DO
        YMIN=V((NREFBA-1)*NLINBL+1)/FSCALE
        YMAX=V(NREFBA*NLINBL)/FSCALE
        DO L=1,NLINBA
          DO I=1,NXYMAX
            SX=U((L-1)*NLINBL+NREFBL)/FSCALE
            SY=YMIN+REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)
            CALL FMAP(N,AIJ_,BIJ_,SX,SY,XP(I),YP(I))
          END DO
          CALL PGSCI(5)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END DO
        CALL PGEBUF
!------------------------------------------------------------------------------
        LMAPPING=.TRUE.
!------------------------------------------------------------------------------
! limites del grid
        SXMINGRID=NINT(XMIN)
        SYMINGRID=NINT(YMIN)
        SXMAXGRID=NINT(XMAX)
        SYMAXGRID=NINT(YMAX)
        WRITE(*,100) '=> SXMINGRID,SYMINGRID: '
        WRITE(*,*) SXMINGRID,SYMINGRID
        WRITE(*,100) '=> SXMAXGRID,SYMAXGRID: '
        WRITE(*,*) SXMAXGRID,SYMAXGRID
        NDELTAX=SXMAXGRID-SXMINGRID
        NDELTAY=SYMAXGRID-SYMINGRID
        WRITE(*,100) '=> Delta X, Delta Y: '
        WRITE(*,*) NDELTAX,NDELTAY
        IF(NDELTAX.GT.NXMAX)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> Delta X > NXMAX'
          INCLUDE 'deallocate_arrays.inc'
          STOP
        END IF
        IF(NDELTAY.GT.NYMAX)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> Delta Y > NYMAX'
          INCLUDE 'deallocate_arrays.inc'
          STOP
        END IF
!------------------------------------------------------------------------------
! si se desea, podemos extrapolar el grid
        CEXT(1:1)=READC('Are you extrapolating the grid (y/n)','n','yn')
        IF(CEXT.EQ.'y')THEN
          SXMINEXTG=READILIM('No. of pixels (X<Xmin)','0',0,NXMAX-NDELTAX)
          SXMAXEXTG=READILIM('No. of pixels (X>Xmax)','0',0,NXMAX-NDELTAX-SXMINEXTG)
          SYMINEXTG=READILIM('No. of pixels (Y<Ymin)','0',0,NYMAX-NDELTAY)
          SYMAXEXTG=READILIM('No. of pixels (Y>Ymax)','0',0,NYMAX-NDELTAY-SYMINEXTG)
        ELSE
          SXMINEXTG=0
          SXMAXEXTG=0
          SYMINEXTG=0
          SYMAXEXTG=0
        END IF
!------------------------------------------------------------------------------
! dibujamos una malla
        CALL PMAPP
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
