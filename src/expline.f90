! Usando una posicion de referencia, esta subrutina busca las lineas
! espectrales.
!!!     SUBROUTINE EXPLINE(XC,YC,COEFFBL00,COEFFBA00)
        SUBROUTINE EXPLINE(XC,COEFFBL00,COEFFBA00)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
! subroutine arguments
!!!     REAL XC,YC
        REAL XC
        REAL COEFFBL00(20),COEFFBA00(20)
!
        INCLUDE 'largest.inc'
        INCLUDE 'dimensions.inc'
!
        INTEGER NMAXPEAKS
        PARAMETER (NMAXPEAKS=NXYMAX/4)        !supongo que esto sera suficiente
!
        INTEGER READI
        REAL READF
        REAL FPOLY
        CHARACTER*255 READC
!
        INTEGER NCBUFF
        INTEGER I,NS,NP,NS_,NP_,NPMIN_(NXYMAX)
        INTEGER L,LL
        INTEGER NDEGBL(NXYMAX),NDEGBA(NXYMAX)  !polynomial degrees for boundary
        INTEGER NDEGBL00,NDEGBA00                 !maximum of NDEGBL and NDEGBA
        INTEGER NLINBL,NLINBA            !number of BL and BA lines in boundary
        INTEGER NLINBA_                     !number of BA lines in new boundary
        INTEGER NYFINDAL,NINLINE,NMED
        INTEGER NPEAKS
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NX1,NX2,NY1,NY2
        INTEGER NBUFPEAK(NMAXPEAKS)
        INTEGER NWIDTH
        REAL COEFFBL00_(20),COEFFBA00_(20)
        REAL XP(NXYMAX),YP(NXYMAX),ZP(NXYMAX)
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XXMIN(1),XXMAX(1),YYMIN(1),YYMAX(1)
        REAL XMINBL(NXYMAX),XMAXBL(NXYMAX)
        REAL YMINBL(NXYMAX),YMAXBL(NXYMAX)
        REAL XMINBA(NXYMAX),XMAXBA(NXYMAX)
        REAL YMINBA(NXYMAX),YMAXBA(NXYMAX)
        REAL COEFFBL(20,NXYMAX),COEFFBA(20,NXYMAX)      !pol. coef. in boundary
        REAL COEFFBA_(20,NXYMAX)                    !pol. coef. in new boundary
        REAL S
        REAL XGRID,YGRID
        REAL XEXPECT,YEXPECT
        REAL XPEAKS(NXYMAX),YPEAKS(NXYMAX),ZPEAKS(NXYMAX)
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL XBUFPEAK(NMAXPEAKS,NXYMAX),YBUFPEAK(NMAXPEAKS,NXYMAX)
        REAL X0,Y0,DMIN,DIST
        REAL XX0(1),YY0(1)
        REAL X0_(NXYMAX)
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL CHISQR
        REAL FMEAN,FSIGMA,FSIGMA0
        CHARACTER*1 CUPDATE,CCONT,CADD
        CHARACTER*50 CDUMMY
        LOGICAL IFPEAKNEW(NMAXPEAKS,NXYMAX)
        LOGICAL LDEBUGLOCAL
        LOGICAL LINSIDE
        LOGICAL LUSE(NXYMAX)           !indica si la nueva linea se acepta o no
!
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKIMAGEN2/NCBUFF
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKBOUND1/NLINBL,NDEGBL,NDEGBL00
        COMMON/BLKBOUND1B/NLINBA,NDEGBA,NDEGBA00
        COMMON/BLKBOUND2/COEFFBL
        COMMON/BLKBOUND2B/COEFFBA
        COMMON/BLKBOUND4/XMINBL,YMINBL,XMAXBL,YMAXBL
        COMMON/BLKBOUND5/XMINBA,YMINBA,XMAXBA,YMAXBA
        COMMON/BLKDEFAULTS1/NYFINDAL,NINLINE
        COMMON/BLKDEFAULTS2/NWIDTH
!------------------------------------------------------------------------------
        LDEBUGLOCAL=.TRUE.
! dibujamos linea en la direccion spectral
        IF(LDEBUGLOCAL)THEN
          XMIN=XMINBL(1)
          XMAX=XMAXBL(1)
          DO I=1,NXYMAX
            XP(I)=XMIN+REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)
            YP(I)=FPOLY(NDEGBL00,COEFFBL00,XP(I))
          END DO
          CALL PGSCI(7)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END IF
! intersecciones de la linea local con los espectros del borde
        CALL INTERSEC(NDEGBA00,COEFFBA00,NDEGBL(1),COEFFBL(1,1),XC,XMIN,YMIN)
        CALL INTERSEC(NDEGBA00,COEFFBA00,NDEGBL(NLINBL),COEFFBL(1,NLINBL),XC,XMAX,YMAX)
! dibujamos linea en la direccion espacial
        IF(LDEBUGLOCAL)THEN
          DO I=1,NXYMAX
            YP(I)=YMIN+REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)
            XP(I)=FPOLY(NDEGBA00,COEFFBA00,YP(I))
          END DO
          CALL PGSCI(7)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END IF
!------------------------------------------------------------------------------
! Grid en Y con NYFINDAL lineas: buscamos los picos en los NYFINDAL espectros
        DO NS=1,NYFINDAL
          S=REAL(NS)/REAL(NYFINDAL+1)
          CALL X2ARC(NDEGBA00,COEFFBA00,YMIN,YMAX,S,YGRID)
          XGRID=FPOLY(NDEGBA00,COEFFBA00,YGRID)
          CALL WHEREAMI(XGRID,YGRID,COEFFBL00_,COEFFBA00_,3,-1,.FALSE.)
          DO I=1,NAXIS(1,NCBUFF)
            XP(I)=REAL(I)
            YP(I)=FPOLY(NDEGBL00,COEFFBL00_,XP(I))
            ZP(I)=IMAGEN(I,NINT(YP(I)),NCBUFF)
          END DO
          CALL FINDMAX(NAXIS(1,NCBUFF),1,NAXIS(1,NCBUFF),XP,ZP,NPEAKS,XPEAKS,ZPEAKS)
          IF(NPEAKS.GT.NMAXPEAKS)THEN
            WRITE(*,101) '***FATAL ERROR***'
            WRITE(*,100) '=> NPEAKS, NMAXPEAKS: '
            WRITE(*,*) NPEAKS,NMAXPEAKS
            WRITE(*,101) '=> NPEAKS > NMAXPEAKS in EXTLINE.'
            WRITE(*,101) '=> Redefine NMAXPEAKS parameter.'
            INCLUDE 'deallocate_arrays.inc'
            STOP
          END IF
          IF(LDEBUGLOCAL)THEN
            CALL SUBPLOT(NAXIS(1,NCBUFF),NX1,NX2,XP,ZP,XP,ZP,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
             'wavelength direction','signal',' ',3,201,1.0)
            CALL SUBPLOTBIS(NPEAKS,1,NPEAKS,XPEAKS,ZPEAKS,XPEAKS,ZPEAKS,.FALSE.,.FALSE.,2,21,1.5)
          END IF
          DO I=1,NPEAKS
            YPEAKS(I)=FPOLY(NDEGBL00,COEFFBL00_,XPEAKS(I))
          END DO
          CALL PGSCI(2)
          CALL PGPOINT(NPEAKS,XPEAKS,YPEAKS,21)
          CALL PGSCI(1)
          NBUFPEAK(NS)=NPEAKS                !numero de picos en el espectro NS
          DO I=1,NPEAKS                   !almacenamos la posicion de los picos
            XBUFPEAK(I,NS)=XPEAKS(I)
            YBUFPEAK(I,NS)=YPEAKS(I)
            IFPEAKNEW(I,NS)=.TRUE.
          END DO
        END DO
        CCONT(1:1)=READC('Do you want to continue (y/n)','y','yn')
        IF(CCONT.EQ.'n') RETURN
!------------------------------------------------------------------------------
! Reconocemos las lineas. Para ello, vamos recorriendo los NYFINDAL espectros,
! comprobando si un nuevo pico (no utilizado anteriormente) tiene nuevas
! identificaciones en el resto de los NYFINDAL-1 espectros.
!------------------------------------------------------------------------------
! Dado que exigimos tener al menos NINLINE picos detectados, no tenemos que
! hace el primer bucle que sigue a continuacion desde 1 hasta NYFINDAL, sino 
! que basta hacerlo desde 1 hasta NYFINDAL-NINLINE+1. Llegados al espectro
! numero NYFINDAL-NINLINE+2, quedan NINLINE-1 espectros por examinar, con lo
! cual nunca podremos obtener NINLINE picos.
!
        NLINBA_=0  !recalculamos las lineas espectrales que definen el boundary
        DO NS=1,NYFINDAL-NINLINE+1  !recorremos los espectros de abajo a arriba
          DO NP=1,NBUFPEAK(NS)     !recorremos los picos de izquierda a derecha
! tomamos el siguiente pico no utilizado todavia
            IF(IFPEAKNEW(NP,NS))THEN
              IFPEAKNEW(NP,NS)=.FALSE.   !ha dejado de ser un pico sin utilizar
              X0=XBUFPEAK(NP,NS)
              Y0=YBUFPEAK(NP,NS)
              IF(LDEBUGLOCAL)THEN
                !usamos un array unidimensional porque el compilador
                !gfortran-mp-10 da error al usar un escalar en lugar
                !de una matriz
                XX0(1)=X0
                YY0(1)=Y0
                CALL PGSCI(5)
                CALL PGPOINT(1,XX0,YY0,17)
                CALL PGSCI(1)
              END IF
! determinamos condiciones locales
              CALL WHEREAMI(X0,Y0,COEFFBL00_,COEFFBA00_,-1,-1,.FALSE.)
! calculamos extremos de la linea espectral
              CALL INTERSEC(NDEGBA00,COEFFBA00_,NDEGBL(1),COEFFBL(1,1),X0,XMIN,YMIN)
              CALL INTERSEC(NDEGBA00,COEFFBA00_,NDEGBL(NLINBL),COEFFBL(1,NLINBL),X0,XMAX,YMAX)
              IF(LDEBUGLOCAL)THEN
                !usamos un array unidimensional porque el compilador
                !gfortran-mp-10 da error al usar un escalar en lugar
                !de una matriz
                XXMIN(1)=XMIN
                XXMAX(1)=XMAX
                YYMIN(1)=YMIN
                YYMAX(1)=YMAX
                CALL PGSCI(7)
                CALL PGPOINT(1,XXMIN,YYMIN,17)
                CALL PGPOINT(1,XXMAX,YYMAX,17)
                CALL PGSCI(1)
              END IF
! chequeamos si los limites estan dentro de la imagen (con un borde de
! seguridad determinado por NMED)
              NMED=3
              LINSIDE=((XMIN.GT.REAL(1+NMED)).AND.(XMIN.LT.REAL(NAXIS(1,NCBUFF)-NMED)).AND. &
                       (YMIN.GT.REAL(1+NMED)).AND.(YMIN.LT.REAL(NAXIS(2,NCBUFF)-NMED)).AND. &
                       (XMAX.GT.REAL(1+NMED)).AND.(XMAX.LT.REAL(NAXIS(1,NCBUFF)-NMED)).AND. &
                       (YMAX.GT.REAL(1+NMED)).AND.(YMAX.LT.REAL(NAXIS(2,NCBUFF)-NMED)))
! OK! estamos dentro
              IF(LINSIDE)THEN
!..............................................................................
! Predecimos posiciones esperadas de la linea en los demas NYFINDLA espectros,
! y comprobamos si hay algun pico cerca; la linea es aceptada si el numero
! de picos proximos es igual o mayor que NINLINE. Como estamos barriendo los
! espectros de abajo a arriba, no es necesario chequear los picos de los
! espectros que quedan por debajo del actual (esos picos ya han sido examinados
! y si no han sido ubicados en ninguna linea, tampoco deben ser utiles para
! los picos que se encuentran en el espectro actual y por encima). 
                NPEAKS=1     !tenemos garantizado el pico en el espectro actual
                XF(NPEAKS)=X0
                YF(NPEAKS)=Y0
                DO NS_=NS+1,NYFINDAL
                  S=REAL(NS_)/REAL(NYFINDAL+1)
                  CALL X2ARC(NDEGBA00,COEFFBA00_,YMIN,YMAX,S,YEXPECT)
                  XEXPECT=FPOLY(NDEGBA00,COEFFBA00_,YEXPECT)
                  NPMIN_(NS_)=0           !vamos a buscar el pico mas cercano
                  DMIN=10.*REAL(NXYMAX)
                  DO NP_=1,NBUFPEAK(NS_)
                    IF(IFPEAKNEW(NP_,NS_))THEN             !pico sin usar aun
                      DIST=ABS(XBUFPEAK(NP_,NS_)-XEXPECT)
                      IF(DIST.LT.DMIN)THEN
                        DMIN=DIST
                        NPMIN_(NS_)=NP_
                      END IF
                    END IF
                  END DO
                  IF(NPMIN_(NS_).GT.0)THEN       !quedaba algun pico sin usar
                    IF(NINT(DMIN).LE.NMED)THEN  !si esta cerca de lo esperado
                      NPEAKS=NPEAKS+1
                      XF(NPEAKS)=XBUFPEAK(NPMIN_(NS_),NS_)
                      YF(NPEAKS)=YBUFPEAK(NPMIN_(NS_),NS_)
                    END IF
                  END IF
                END DO
! si el numero de picos encontrados es suficientemente grande, ajustamos la
! linea y eliminamos los picos usados del pool de picos sin usar
                IF(NPEAKS.GE.NINLINE)THEN
                  NLINBA_=NLINBA_+1
                  CALL POLFIT(YF,XF,XF,NPEAKS,NDEGBA00+1,0,COEFFBA_(1,NLINBA_),CHISQR)
                  IF(LDEBUGLOCAL)THEN
                    CALL INTERSEC(NDEGBA00,COEFFBA_(1,NLINBA_),NDEGBL(1),COEFFBL(1,1),X0,XMIN,YMIN)
                    CALL INTERSEC(NDEGBA00,COEFFBA_(1,NLINBA_),NDEGBL(NLINBL),COEFFBL(1,NLINBL),X0,XMAX,YMAX)
                    DO I=1,NXYMAX
                      YP(I)=YMIN+REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)
                      XP(I)=FPOLY(NDEGBA00,COEFFBA_(1,NLINBA_),YP(I))
                    END DO
                    CALL PGSCI(5)
                    CALL PGLINE(NXYMAX,XP,YP)
                    CALL PGSCI(1)
                    WRITE(CDUMMY,*) NLINBA_
                    CALL RMBLANK(CDUMMY,CDUMMY,L)
                    CALL SUBPLOT(NPEAKS,1,NPEAKS,YF,XF,XF,YF,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
                     'Y axis','X axis','Line #'//CDUMMY(1:L),2,17,1.5)
                    CALL SUBPLOTBIS(NXYMAX,1,NXYMAX,YP,XP,YP,XP,.FALSE.,.FALSE.,3,101,1.0)
                  END IF
                  X0_(NLINBA_)=X0           !hace falta para refinar las lineas
                  DO NS_=NS,NYFINDAL !eliminamos los picos del pool
                    IF(NS_.EQ.NS)THEN
                      IFPEAKNEW(NP,NS_)=.FALSE.
                    ELSE
                      IF(NPMIN_(NS_).NE.0)THEN
                        IFPEAKNEW(NPMIN_(NS_),NS_)=.FALSE.
                      END IF
                    END IF
                  END DO
                END IF
!..............................................................................
              END IF !comprueba si la posible linea esta entera en la IMAGEN
            END IF !comprueba si IFPEAKNEW(NP,NS) es .TRUE.
          END DO !bucle en numero de picos NP
        END DO !buble en numero de espectros NS
!------------------------------------------------------------------------------
        WRITE(*,101) '***WARNING***'
        WRITE(*,100) '=> previous number of BA lines: '
        WRITE(*,*) NLINBA
        WRITE(*,100) '=> initial number of new lines: '
        WRITE(*,*) NLINBA_
        WRITE(CDUMMY,*) NWIDTH
        NWIDTH=READI('No. of pixels to fit gaussians',CDUMMY)
        I=0
        FSIGMA0=0.07
        DO L=1,NLINBA_
          WRITE(CDUMMY,'(I10,A1,I10)') L,',',NLINBA_
          CALL RMBLANK(CDUMMY,CDUMMY,LL)
          WRITE(*,100) '=> Refining line #'
          WRITE(*,101) CDUMMY(1:LL)
          CALL REFINEBA(NDEGBA00,COEFFBA_(1,L),X0_(L),FMEAN,FSIGMA)
          IF(CCONT.EQ.'a')THEN
            LUSE(L)=(FSIGMA.LE.FSIGMA0)
          ELSE
            IF(FSIGMA.LE.FSIGMA0)THEN
              CCONT='y'
            ELSE
              CCONT='n'
            END IF
            CCONT(1:1)=READC('Are you using this last fit (y/n/a=auto)',CCONT,'yna')
            IF(CCONT.EQ.'a')THEN
              WRITE(CDUMMY,*) FSIGMA0
              FSIGMA0=READF('Sigma threshold',CDUMMY)
              LUSE(L)=(FSIGMA.LE.FSIGMA0)
            ELSE
              LUSE(L)=(CCONT.EQ.'y')
            END IF
          END IF
          IF(LUSE(L)) I=I+1
        END DO
        WRITE(*,100) '=> final number of new lines..: '
        WRITE(*,*) I
        IF(I.EQ.0) RETURN
!------------------------------------------------------------------------------
! si nos gusta el resultado, podemos actualizar las lineas BA que definen
! el boundary
        CUPDATE(1:1)=READC('Do you want to update BA lines (y/n)','y','yn')
        IF(CUPDATE.EQ.'y')THEN
          CADD(1:1)=READC('[a]ppend to previous list or define [n]ew list','a','an')
          IF(CADD.EQ.'n') NLINBA=0
          DO L=1,NLINBA_
            IF(LUSE(L))THEN
              NLINBA=NLINBA+1
              NDEGBA(NLINBA)=NDEGBA00
              DO I=1,NDEGBA00+1
                COEFFBA(I,NLINBA)=COEFFBA_(I,L)
              END DO
            END IF
          END DO
          CALL SORTLINES !esto es fundamental, las lineas no estan ordenadas
          WRITE(*,101) '***WARNING***'
          WRITE(*,101) '=> RESET shrinking/stretching'
          WRITE(*,*)
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
