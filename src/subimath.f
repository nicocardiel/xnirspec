C Version 06-July-2007
C calculadora
C
        SUBROUTINE SUBIMATH
        IMPLICIT NONE
C
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
C
        INTEGER READILIM
        INTEGER TRUELEN
        INTEGER SYSTEMFUNCTION
        REAL READF
        REAL FMEAN0,FMEAN0E,FMEAN1,FMEAN2,FMEAN2E,FMEDIAN1
        REAL FROBUSTRMS0,FROBUSTRMS1
        CHARACTER*255 READC
C
        INTEGER I,J,K,L,I1,I2,J1,J2,I0,J0
        INTEGER II,JJ
        INTEGER NORIGEN,NB,NB_
        INTEGER NB_CANCEL,NB_EXIT
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NFRAMES(NMAXBUFF)
        INTEGER NMED1,NMED2
        INTEGER NQUAD0,NQUAD1
        INTEGER NOFFX0,NOFFX1,NOFFY0,NOFFY1
        INTEGER NBUFF0,NBUFF1
        INTEGER NX1,NX2,NY1,NY2
        INTEGER NX1_,NX2_,NY1_,NY2_
        INTEGER ICUT,IWIDTH
        INTEGER NBOX_X,NBOX_Y
        INTEGER NEXTINFO
        INTEGER NDEGREE,NCOEF,NCOEFDER
        INTEGER DI(9),DJ(9)
        INTEGER ISYSTEM
        INTEGER NX1_PLOT,NX2_PLOT,NY1_PLOT,NY2_PLOT
        INTEGER IKERNEL
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        REAL PIXEL(NXYMAX*NXYMAX) !OJO: para no consumir mas memoria, esta
                                  !     variable utiliza el mismo common que
                                  !     la variable IMAGEN_
        REAL XCUT(NXMAX),YCUT(NYMAX)
        REAL XCUT_ERR(NXMAX),YCUT_ERR(NYMAX)
        REAL XC,YC
        REAL FCTE,FCTE_ERR,FZERO,FZERO_ERR
        REAL FSIGMA,TSIGMA
        REAL XF(NXMAX),YF(NYMAX)
        REAL COEF(20),CHISQR
        REAL COEFDER(20)
        REAL FACTOR,SUMFACTOR
        REAL XOFFSET,YOFFSET,FI0,FJ0
        REAL THRESHOLD
        REAL NEWCONSTANT,THRESHOLD1,THRESHOLD2
        REAL FDUMMY
        REAL OFFSET_LOG
        CHARACTER*1 CH,CZERO
        CHARACTER*1 COPER
        CHARACTER*1 CUTIL
        CHARACTER*1 CFILT,CAXIS,CFUN_OR_DER,CRMS
        CHARACTER*1 CKERNEL
        CHARACTER*50 CDUMMY
        LOGICAL LWAVECAL(NMAXBUFF)
        LOGICAL LDEFBUFF(NMAXBUFF)
        LOGICAL LEXIT,LCANCEL,LCONT,LASK
        LOGICAL LDIVZERO
        LOGICAL IFCHAN(NXMAX),IFSCAN(NYMAX)
        LOGICAL LFIT(NXYMAX)
        LOGICAL LOK
C
        COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKIMAGEN1_/PIXEL                !es global para ahorrar memoria
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1         !wavelength calibration
        COMMON/BLKWAVECAL2/LWAVECAL                     !wavelength calibration
        COMMON/BLKNAXIS/NAXIS                                      !dimensiones
        COMMON/BLKNFRAMES/NFRAMES
        COMMON/BLKDEFAULTS3/TSIGMA
        COMMON/BLKXYLIMPLOT/NX1_PLOT,NX2_PLOT,NY1_PLOT,NY2_PLOT
C------------------------------------------------------------------------------
C Note: the pattern of the frames in box-9 is the following:
C       6 9 4
C       3 1 7
C       8 5 2
C offsets of each frame in the 768x768 composite mosaic
        DATA (DI(K),DJ(K),K=1,9) /
     +   256,256,  !1
     +   000,512,  !2
     +   256,000,  !3
     +   512,512,  !4
     +   000,256,  !5
     +   512,000,  !6
     +   256,512,  !7
     +   000,000,  !8
     +   512,256/  !9
C
        I2=0 !evita WARNING de compilacion
        J2=0 !evita WARNING de compilacion
        NBUFF0=0 !evita WARNING de compilacion
        NBUFF1=0 !evita WARNING de compilacion
        FCTE_ERR=0.0 !evita WARNING de compilacion
        THRESHOLD=0.0 !evita WARNING de compilacion
        NEWCONSTANT=0.0 !evita WARNING de compilacion
        THRESHOLD1=0.0 !evita WARNING de compilacion
        THRESHOLD2=0.0 !evita WARNING de compilacion
C------------------------------------------------------------------------------
        NORIGEN=86
        CFILT='@'
        CFUN_OR_DER='@'
        CRMS='@'
        IWIDTH=3
        NBOX_X=3
        NBOX_Y=3
        NDEGREE=0
        CAXIS='@'
        LEXIT=.FALSE.
        DO WHILE(.NOT.LEXIT)
          CALL RPGERASW(0.49,0.91,0.10,0.70,4)
          LCANCEL=.FALSE.
C..............................................................................
C elegimos primer buffer
          DO I=1,NMAXBUFF
            NB=10*(I-1)+NORIGEN
            IF(I.LE.NMAXBUFF/2)THEN
              WRITE(CDUMMY,'(A2,I1,A1)') '#[',I,']'
            ELSE
              WRITE(CDUMMY,'(A5,I1)') 'err #',I-NMAXBUFF/2
            END IF
            L=TRUELEN(CDUMMY)
            CALL BUTTON(NB,CDUMMY(1:L),0)
            LDEFBUFF(I)=((NAXIS(1,I).NE.0).AND.(NAXIS(2,I).NE.0))
            IF(LDEFBUFF(I)) CALL BUTTON(NB,CDUMMY(1:L),1)
          END DO
          NB_CANCEL=10*(NMAXBUFF+1)+NORIGEN
          NB_EXIT=NB_CANCEL+1
          CALL BUTTON(NB_CANCEL,'[C]ancel',0)
          CALL BUTTON(NB_EXIT,'e[X]it',0)
          LCONT=.FALSE.
          DO WHILE(.NOT.LCONT)
            CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
            CALL IFBUTTON(XC,YC,NB)
            IF(CH.EQ.'C')THEN
              WRITE(*,101) '[CANCEL]'
              WRITE(*,*)
              NB=NB_CANCEL
            ELSEIF(CH.EQ.'X')THEN
              WRITE(*,101) '[EXIT]'
              WRITE(*,*)
              NB=NB_EXIT
            ELSEIF(INDEX('123456',CH).NE.0)THEN
              READ(CH,*) NB_
              NB=10*(NB_-1)+NORIGEN
            END IF
            IF(NB.EQ.NB_EXIT)THEN
              LEXIT=.TRUE.
              LCANCEL=.TRUE.
              LCONT=.TRUE.
            ELSEIF(NB.EQ.NB_CANCEL)THEN
              LCANCEL=.TRUE.
              LCONT=.TRUE.
            ELSE
              DO I=1,NMAXBUFF
                IF(NB.EQ.NORIGEN+(I-1)*10)THEN
                  NBUFF0=I
                  IF(I.LE.NMAXBUFF/2)THEN
                    WRITE(CDUMMY,'(A2,I1,A1)') '#[',I,']'
                    WRITE(*,*)
                    WRITE(*,100) 'Selecting buffer #['
                    WRITE(*,'(I1,$)') I
                    WRITE(*,101) ']'
                  ELSE
                    WRITE(CDUMMY,'(A5,I1)') 'err #',I-NMAXBUFF/2
                  END IF
                  L=TRUELEN(CDUMMY)
                  CALL BUTTON(NB,CDUMMY(1:L),5)
                  CALL BUTTON(NB,CDUMMY(1:L),-4)
                  LCONT=.TRUE.
                END IF
              END DO
            END IF
          END DO
          IF(.NOT.LCANCEL)THEN
            DO I=1,NMAXBUFF
              IF(I.NE.NBUFF0)THEN
                IF(I.LE.NMAXBUFF/2)THEN
                  WRITE(CDUMMY,'(A2,I1,A1)') '#[',I,']'
                ELSE
                  WRITE(CDUMMY,'(A5,I1)') 'err #',I-NMAXBUFF/2
                END IF
                L=TRUELEN(CDUMMY)
                NB=10*(I-1)+NORIGEN
                CALL BUTTON(NB,CDUMMY(1:L),0)
                CALL BUTTON(NB,CDUMMY(1:L),3)
              END IF
            END DO
          END IF
C..............................................................................
C mostramos posibles operaciones
          COPER=' '  !evita un WARNING de compilacion (y un posible error)
          IF(.NOT.LCANCEL)THEN
            CALL BUTTON(NORIGEN+01,'# [<] #',0)
            CALL BUTTON(NORIGEN+11,'[+]',0)
            CALL BUTTON(NORIGEN+21,'[-]',0)
            CALL BUTTON(NORIGEN+31,'[*]',0)
            CALL BUTTON(NORIGEN+41,'[/]',0)
            CALL BUTTON(NORIGEN+51,'[0]',0)
            CALL BUTTON(NORIGEN+61,'[f]ilter',0)
            CALL BUTTON(NORIGEN+71,'< 1 -> [4]',0)
            CALL BUTTON(NORIGEN+81,'< 4 -> [1]',0)
            CALL BUTTON(NORIGEN+91,'< [q]uad.',0)
            CALL BUTTON(NORIGEN+101,'[s]hift xy',0)
            CALL BUTTON(NORIGEN+111,'< 1 -> 9',0)
            CALL BUTTON(NORIGEN+121,'< 9 -> 1',0)
            LCONT=.FALSE.
            DO WHILE(.NOT.LCONT)
              CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
              CALL IFBUTTON(XC,YC,NB)
              IF(CH.EQ.'C')THEN
                NB=NB_CANCEL
              ELSEIF(CH.EQ.'X')THEN
                NB=NB_EXIT
              ELSEIF(CH.EQ.'<')THEN
                NB=NORIGEN+01
              ELSEIF(CH.EQ.'+')THEN
                NB=NORIGEN+11
              ELSEIF(CH.EQ.'-')THEN
                NB=NORIGEN+21
              ELSEIF(CH.EQ.'*')THEN
                NB=NORIGEN+31
              ELSEIF(CH.EQ.'/')THEN
                NB=NORIGEN+41
              ELSEIF(CH.EQ.'0')THEN
                NB=NORIGEN+51
              ELSEIF(CH.EQ.'f')THEN
                NB=NORIGEN+61
              ELSEIF(CH.EQ.'4')THEN
                NB=NORIGEN+71
              ELSEIF(CH.EQ.'1')THEN
                NB=NORIGEN+81
              ELSEIF(CH.EQ.'q')THEN
                NB=NORIGEN+91
              ELSEIF(CH.EQ.'s')THEN
                NB=NORIGEN+101
              END IF
              LCONT=.TRUE.
              IF(NB.EQ.NB_EXIT)THEN
                WRITE(*,101) '[EXIT]'
                WRITE(*,*)
                LEXIT=.TRUE.
                LCANCEL=.TRUE.
              ELSEIF(NB.EQ.NB_CANCEL)THEN
                WRITE(*,101) '[CANCEL]'
                WRITE(*,*)
                LCANCEL=.TRUE.
              ELSEIF(NB.EQ.NORIGEN+01)THEN
                COPER='<'
                WRITE(*,100) '['
                WRITE(*,100) COPER
                WRITE(*,101) ']'
              ELSEIF(NB.EQ.NORIGEN+11)THEN
                COPER='+'
                WRITE(*,100) '['
                WRITE(*,100) COPER
                WRITE(*,101) ']'
              ELSEIF(NB.EQ.NORIGEN+21)THEN
                COPER='-'
                WRITE(*,100) '['
                WRITE(*,100) COPER
                WRITE(*,101) ']'
              ELSEIF(NB.EQ.NORIGEN+31)THEN
                COPER='*'
                WRITE(*,100) '['
                WRITE(*,100) COPER
                WRITE(*,101) ']'
              ELSEIF(NB.EQ.NORIGEN+41)THEN
                COPER='/'
                WRITE(*,100) '['
                WRITE(*,100) COPER
                WRITE(*,101) ']'
              ELSEIF(NB.EQ.NORIGEN+51)THEN
                COPER='0'
                WRITE(*,100) '['
                WRITE(*,100) COPER
                WRITE(*,101) ']'
              ELSEIF(NB.EQ.NORIGEN+61)THEN
                COPER='f'
                WRITE(*,101) '[f]ilter'
              ELSEIF(NB.EQ.NORIGEN+71)THEN
                COPER='4'
                WRITE(*,101) '< 1 -> [4]'
              ELSEIF(NB.EQ.NORIGEN+81)THEN
                COPER='1'
                WRITE(*,101) '< 4 -> [1]'
              ELSEIF(NB.EQ.NORIGEN+91)THEN
                COPER='q'
                WRITE(*,101) '< [q]uad.'
              ELSEIF(NB.EQ.NORIGEN+101)THEN
                COPER='s'
                WRITE(*,101) '[s]hift xy'
              ELSEIF(NB.EQ.NORIGEN+111)THEN
                COPER='9'
                WRITE(*,101) '< 1 -> 9'
              ELSEIF(NB.EQ.NORIGEN+121)THEN
                COPER='c'
                WRITE(*,101) '< 9 -> 1'
              ELSE
                LCONT=.FALSE.
              END IF
            END DO
            CALL BUTTON(NORIGEN+01,'# [<] #',3)
            CALL BUTTON(NORIGEN+11,'[+]',3)
            CALL BUTTON(NORIGEN+21,'[-]',3)
            CALL BUTTON(NORIGEN+31,'[*]',3)
            CALL BUTTON(NORIGEN+41,'[/]',3)
            CALL BUTTON(NORIGEN+51,'[0]',3)
            CALL BUTTON(NORIGEN+61,'[f]ilter',3)
            CALL BUTTON(NORIGEN+71,'< 1 -> [4]',3)
            CALL BUTTON(NORIGEN+81,'< 4 -> [1]',3)
            CALL BUTTON(NORIGEN+91,'< [q]uad.',3)
            CALL BUTTON(NORIGEN+101,'[s]hift xy',3)
            CALL BUTTON(NORIGEN+111,'< 1 -> 9',3)
            CALL BUTTON(NORIGEN+121,'< 9 -> 1',3)
            IF(COPER.EQ.'<')THEN
              CALL BUTTON(NORIGEN+01,'# [<] #',5)
              CALL BUTTON(NORIGEN+01,'# [<] #',-4)
            ELSEIF(COPER.EQ.'+')THEN
              CALL BUTTON(NORIGEN+11,'[+]',5)
              CALL BUTTON(NORIGEN+11,'[+]',-4)
            ELSEIF(COPER.EQ.'-')THEN
              CALL BUTTON(NORIGEN+21,'[-]',5)
              CALL BUTTON(NORIGEN+21,'[-]',-4)
            ELSEIF(COPER.EQ.'*')THEN
              CALL BUTTON(NORIGEN+31,'[*]',5)
              CALL BUTTON(NORIGEN+31,'[*]',-4)
            ELSEIF(COPER.EQ.'/')THEN
              CALL BUTTON(NORIGEN+41,'[/]',5)
              CALL BUTTON(NORIGEN+41,'[/]',-4)
            ELSEIF(COPER.EQ.'0')THEN
              CALL BUTTON(NORIGEN+51,'[0]',5)
              CALL BUTTON(NORIGEN+51,'[0]',-4)
            ELSEIF(COPER.EQ.'f')THEN
              CALL BUTTON(NORIGEN+61,'[f]ilter',5)
              CALL BUTTON(NORIGEN+61,'[f]ilter',-4)
            ELSEIF(COPER.EQ.'4')THEN
              CALL BUTTON(NORIGEN+71,'< 1 -> [4]',5)
              CALL BUTTON(NORIGEN+71,'< 1 -> [4]',-4)
            ELSEIF(COPER.EQ.'1')THEN
              CALL BUTTON(NORIGEN+81,'< 4 -> [1]',5)
              CALL BUTTON(NORIGEN+81,'< 4 -> [1]',-4)
            ELSEIF(COPER.EQ.'q')THEN
              CALL BUTTON(NORIGEN+91,'< [q]uad.',5)
              CALL BUTTON(NORIGEN+91,'< [q]uad.',-4)
            ELSEIF(COPER.EQ.'s')THEN
              CALL BUTTON(NORIGEN+101,'[s]hift xy',5)
              CALL BUTTON(NORIGEN+101,'[s]hift xy',-4)
            ELSEIF(COPER.EQ.'9')THEN
              CALL BUTTON(NORIGEN+111,'< 1 -> 9',5)
              CALL BUTTON(NORIGEN+111,'< 1 -> 9',-4)
            ELSEIF(COPER.EQ.'c')THEN
              CALL BUTTON(NORIGEN+121,'< 9 -> 1',5)
              CALL BUTTON(NORIGEN+121,'< 9 -> 1',-4)
            END IF
          END IF
C..............................................................................
C elegimos region de la imagen inicial a utilizar (salvo si estamos copiando
C un buffer completo o haciendo una transformacion 1 <--> 4 o 1 <--> 9, en 
C cuyo caso no hay nada que elegir)
          IF(.NOT.LCANCEL)THEN
            CALL BUTTON(NORIGEN+02,'[f]ull frame',0)
            CALL BUTTON(NORIGEN+12,'[r]egion',0)
            LCONT=.FALSE.
            IF((COPER.EQ.'<').OR.(COPER.EQ.'4').OR.
     +       (COPER.EQ.'1').OR.(COPER.EQ.'9').OR.
     +       (COPER.EQ.'c').OR.(COPER.EQ.'q'))THEN
              NX1=0
              NX2=0
              NY1=0
              NY2=0
              LCONT=.TRUE.
              CALL BUTTON(NORIGEN+02,'[f]ull frame',5)
              CALL BUTTON(NORIGEN+02,'[f]ull frame',-4)
              CALL BUTTON(NORIGEN+12,'[r]egion',3)
            END IF
            DO WHILE(.NOT.LCONT)
              CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
              CALL IFBUTTON(XC,YC,NB)
              IF(CH.EQ.'C')THEN
                NB=NB_CANCEL
              ELSEIF(CH.EQ.'X')THEN
                NB=NB_EXIT
              ELSEIF(CH.EQ.'f')THEN
                NB=NORIGEN+02
              ELSEIF(CH.EQ.'r')THEN
                NB=NORIGEN+12
              END IF
              LCONT=.TRUE.
              IF(NB.EQ.NB_EXIT)THEN
                WRITE(*,101) '[EXIT]'
                WRITE(*,*)
                LEXIT=.TRUE.
                LCANCEL=.TRUE.
              ELSEIF(NB.EQ.NB_CANCEL)THEN
                WRITE(*,101) '[CANCEL]'
                WRITE(*,*)
                LCANCEL=.TRUE.
              ELSEIF(NB.EQ.NORIGEN+02)THEN
                WRITE(*,101) '[f]ull frame'
                CALL BUTTON(NORIGEN+02,'[f]ull frame',5)
                CALL BUTTON(NORIGEN+02,'[f]ull frame',-4)
                CALL BUTTON(NORIGEN+12,'[r]egion',3)
                NX1=1
                NX2=NAXIS(1,NBUFF0)
                NY1=1
                NY2=NAXIS(2,NBUFF0)
              ELSEIF(NB.EQ.NORIGEN+12)THEN
                WRITE(*,101) '[r]egion'
                CALL BUTTON(NORIGEN+12,'[r]egion',5)
                CALL BUTTON(NORIGEN+12,'[r]egion',-4)
                CALL BUTTON(NORIGEN+02,'[f]ull frame',3)
                WRITE(CDUMMY,*) NX1_PLOT
                NX1=READILIM('X min (pixel)',CDUMMY,1,NAXIS(1,NBUFF0))
                WRITE(CDUMMY,*) NX2_PLOT
                NX2=READILIM('X max (pixel)',CDUMMY,NX1,NAXIS(1,NBUFF0))
                WRITE(CDUMMY,*) NY1_PLOT
                NY1=READILIM('Y min (pixel)',CDUMMY,1,NAXIS(2,NBUFF0))
                WRITE(CDUMMY,*) NY2_PLOT
                NY2=READILIM('Y max (pixel)',CDUMMY,NY1,NAXIS(2,NBUFF0))
              ELSE
                LCONT=.FALSE.
              END IF
            END DO
          END IF
C..............................................................................
C la operacion '0' puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'0')THEN
              DO I=NY1,NY2
                DO J=NX1,NX2
                  IMAGEN(J,I,NBUFF0)=0.
                END DO
              END DO
              IF(NBUFF0.LE.NMAXBUFF/2)THEN
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=0.
                  END DO
                END DO
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C pedimos segundo elemento de la operacion (salvo si la operacion es copiar
C frame, filtrar frame, shift, 1 <--> 4 o s <--> 9, en cuyo caso esta opcion es
C inmediata)
          CUTIL=' ' !Evita un WARNING de compilacion (y un posible error)
          IF(.NOT.LCANCEL)THEN
            CALL BUTTON(NORIGEN+22,'cons[t]ant',0)
            CALL BUTTON(NORIGEN+32,'[f]rame',0)
            CALL BUTTON(NORIGEN+42,'[x]-cut',0)
            CALL BUTTON(NORIGEN+52,'[y]-cut',0)
            LCONT=.FALSE.
            IF(COPER.EQ.'<')THEN
              CUTIL='f'
              LCONT=.TRUE.
            ELSEIF(COPER.EQ.'4')THEN
              CUTIL='f'
              LCONT=.TRUE.
            ELSEIF(COPER.EQ.'1')THEN
              CUTIL='f'
              LCONT=.TRUE.
            ELSEIF(COPER.EQ.'9')THEN
              CUTIL='f'
              LCONT=.TRUE.
            ELSEIF(COPER.EQ.'c')THEN
              CUTIL='f'
              LCONT=.TRUE.
            ELSEIF(COPER.EQ.'q')THEN
              CUTIL='f'
              LCONT=.TRUE.
            ELSEIF(COPER.EQ.'f')THEN
              CUTIL=' '
              LCONT=.TRUE.
            ELSEIF(COPER.EQ.'s')THEN
              CUTIL=' '
              LCONT=.TRUE.
            END IF
            DO WHILE(.NOT.LCONT)
              CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
              CALL IFBUTTON(XC,YC,NB)
              IF(CH.EQ.'C')THEN
                NB=NB_CANCEL
              ELSEIF(CH.EQ.'X')THEN
                NB=NB_EXIT
              ELSEIF(CH.EQ.'t')THEN
                NB=NORIGEN+22
              ELSEIF(CH.EQ.'f')THEN
                NB=NORIGEN+32
              ELSEIF(CH.EQ.'x')THEN
                NB=NORIGEN+42
              ELSEIF(CH.EQ.'y')THEN
                NB=NORIGEN+52
              END IF
              LCONT=.TRUE.
              IF(NB.EQ.NB_EXIT)THEN
                WRITE(*,101) '[EXIT]'
                WRITE(*,*)
                LEXIT=.TRUE.
                LCANCEL=.TRUE.
              ELSEIF(NB.EQ.NB_CANCEL)THEN
                WRITE(*,101) '[CANCEL]'
                WRITE(*,*)
                LCANCEL=.TRUE.
              ELSEIF(NB.EQ.NORIGEN+22)THEN
                CUTIL='t'
                WRITE(*,101) 'cons[t]ant'
              ELSEIF(NB.EQ.NORIGEN+32)THEN
                CUTIL='f'
                WRITE(*,101) '[f]rame'
              ELSEIF(NB.EQ.NORIGEN+42)THEN
                CUTIL='x'
                WRITE(*,101) '[x] cut'
              ELSEIF(NB.EQ.NORIGEN+52)THEN
                CUTIL='y'
                WRITE(*,101) '[y] cut'
              ELSE
                LCONT=.FALSE.
              END IF
            END DO
            IF(.NOT.LCANCEL)THEN
              CALL BUTTON(NORIGEN+22,'cons[t]ant',3)
              CALL BUTTON(NORIGEN+32,'[f]rame',3)
              CALL BUTTON(NORIGEN+42,'[x]-cut',3)
              CALL BUTTON(NORIGEN+52,'[y]-cut',3)
              IF(CUTIL.EQ.'t')THEN
                CALL BUTTON(NORIGEN+22,'cons[t]ant',5)
                CALL BUTTON(NORIGEN+22,'cons[t]ant',-4)
              ELSEIF(CUTIL.EQ.'f')THEN
                CALL BUTTON(NORIGEN+32,'[f]rame',5)
                CALL BUTTON(NORIGEN+32,'[f]rame',-4)
              ELSEIF(CUTIL.EQ.'x')THEN
                CALL BUTTON(NORIGEN+42,'[x]-cut',5)
                CALL BUTTON(NORIGEN+42,'[x]-cut',-4)
              ELSEIF(CUTIL.EQ.'y')THEN
                CALL BUTTON(NORIGEN+52,'[y]-cut',5)
                CALL BUTTON(NORIGEN+52,'[y]-cut',-4)
              END IF
            END IF
          END IF
C..............................................................................
C NOTA: las operaciones con las imagenes de errores deben hacerse antes que 
C las operaciones con las imagenes de datos, para  evitar que al redefinir la
C matriz de datos, esto provoque que perdamos informacion que luego 
C necesitariamos para calcular las matrices de errores.
C..............................................................................
C la operacion con una constante ya puede hacerse
          IF(.NOT.LCANCEL)THEN
            IF(CUTIL.EQ.'t')THEN
              LASK=.TRUE.
              DO WHILE(LASK)
                FCTE=READF('Constant','@')
                LASK=((FCTE.EQ.0.0).AND.(COPER.EQ.'/'))
                IF(LASK)THEN
                  WRITE(*,100) '=> Invalid Constant'
                  WRITE(*,101) '(division by zero!)'
                  WRITE(*,100) '(press <CR> to continue...)'
                  READ(*,*)
                END IF
                IF(NBUFF0.LE.NMAXBUFF/2)THEN
                  FCTE_ERR=READF('Error in constant','0.0')
                END IF
              END DO
              IF(COPER.EQ.'+')THEN
                IF(NBUFF0.LE.NMAXBUFF/2)THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)+
     +                 FCTE_ERR*FCTE_ERR)
                    END DO
                  END DO
                END IF
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)+FCTE
                  END DO
                END DO
              ELSEIF(COPER.EQ.'-')THEN
                IF(NBUFF0.LE.NMAXBUFF/2)THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)+
     +                 FCTE_ERR*FCTE_ERR)
                    END DO
                  END DO
                END IF
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)-FCTE
                  END DO
                END DO
              ELSEIF(COPER.EQ.'*')THEN
                IF(NBUFF0.LE.NMAXBUFF/2)THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                 IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                 FCTE_ERR*FCTE_ERR+
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*FCTE*FCTE)
                    END DO
                  END DO
                END IF
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)*FCTE
                  END DO
                END DO
              ELSEIF(COPER.EQ.'/')THEN
                IF(NBUFF0.LE.NMAXBUFF/2)THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                 IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                 FCTE_ERR*FCTE_ERR+
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*FCTE*FCTE)/
     +                 (FCTE*FCTE)
                    END DO
                  END DO
                END IF
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)/FCTE
                  END DO
                END DO
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C elegimos segundo buffer (OJO: solo dejamos elegir un buffer ya definido y
C que ademas tenga las mismas dimensiones que el buffer inicial ---esta ultima
C condicion no es exigible cuando estemos copiando un buffer ya definido en 
C uno vacio, cuando estamos filtrando una imagen, cuando comprimimos/expandimos
C un box-9 o cuando hagamos shift)
          IF(.NOT.LCANCEL)THEN
            DO I=1,NMAXBUFF
              NB=10*(I-1)+NORIGEN+3
              IF(I.LE.NMAXBUFF/2)THEN
                WRITE(CDUMMY,'(A2,I1,A1)') '#[',I,']'
              ELSE
                WRITE(CDUMMY,'(A5,I1)') 'err #',I-NMAXBUFF/2
              END IF
              L=TRUELEN(CDUMMY)
              CALL BUTTON(NB,CDUMMY(1:L),0)
              IF((COPER.EQ.'<').OR.(COPER.EQ.'f').OR.
     +         (COPER.EQ.'4').OR.(COPER.EQ.'1').OR.
     +         (COPER.EQ.'9').OR.(COPER.EQ.'c').OR.
     +         (COPER.EQ.'s').OR.(COPER.EQ.'q'))THEN
                LDEFBUFF(I)=((NAXIS(1,I).NE.0).AND.
     +                       (NAXIS(2,I).NE.0))
              ELSE
                LDEFBUFF(I)=((NAXIS(1,I).EQ.NAXIS(1,NBUFF0)).AND.
     +                       (NAXIS(2,I).EQ.NAXIS(2,NBUFF0)))
              END IF
              IF((COPER.EQ.'f').OR.(COPER.EQ.'s'))THEN
                IF(I.EQ.NBUFF0)THEN
                  CALL BUTTON(NB,CDUMMY(1:L),3)
                ELSE
                  IF(LDEFBUFF(I)) CALL BUTTON(NB,CDUMMY(1:L),1)
                END IF
              ELSE
                IF(.NOT.LDEFBUFF(I)) CALL BUTTON(NB,CDUMMY(1:L),3)
              END IF
            END DO
            LCONT=.FALSE.
            DO WHILE(.NOT.LCONT)
              CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
              CALL IFBUTTON(XC,YC,NB)
              IF(CH.EQ.'C')THEN
                NB=NB_CANCEL
              ELSEIF(CH.EQ.'X')THEN
                NB=NB_EXIT
              ELSEIF(INDEX('123456',CH).NE.0)THEN
                READ(CH,*) NB_
                IF((COPER.EQ.'f').OR.(COPER.EQ.'s'))THEN
                  NB=10*(NB_-1)+NORIGEN+3
                ELSE
                  IF(LDEFBUFF(NB_)) NB=10*(NB_-1)+NORIGEN+3
                END IF
              END IF
              IF(NB.EQ.NB_EXIT)THEN
                WRITE(*,101) '[EXIT]'
                LEXIT=.TRUE.
                LCANCEL=.TRUE.
                LCONT=.TRUE.
              ELSEIF(NB.EQ.NB_CANCEL)THEN
                WRITE(*,101) '[CANCEL]'
                LCANCEL=.TRUE.
                LCONT=.TRUE.
              ELSE
                DO I=1,NMAXBUFF
                  IF((NB.EQ.NORIGEN+(I-1)*10+03).AND.
     +             (LDEFBUFF(I).OR.(COPER.EQ.'f').OR.
     +              (COPER.EQ.'s')))THEN
                    IF(COPER.EQ.'9')THEN
                      LOK=((NAXIS(1,I).EQ.256).AND.(NAXIS(2,I).EQ.256))
                    ELSEIF(COPER.EQ.'c')THEN
                      LOK=((NAXIS(1,I).EQ.768).AND.(NAXIS(2,I).EQ.768))
                    ELSE
                      LOK=.TRUE.
                    END IF
                    IF(LOK)THEN
                      NBUFF1=I
                      IF(I.LE.NMAXBUFF/2)THEN
                        WRITE(CDUMMY,'(A2,I1,A1)') '#[',I,']'
                        WRITE(*,100) 'Selecting buffer #['
                        WRITE(*,'(I1,$)') I
                        WRITE(*,101) ']'
                      ELSE
                        WRITE(CDUMMY,'(A5,I1)') 'err #',I-NMAXBUFF/2
                      END IF
                      L=TRUELEN(CDUMMY)
                      CALL BUTTON(NB,CDUMMY(1:L),5)
                      CALL BUTTON(NB,CDUMMY(1:L),-4)
                      LCONT=.TRUE.
                    END IF
                  END IF
                END DO
              END IF
            END DO
          END IF
          IF(.NOT.LCANCEL)THEN
            DO I=1,NMAXBUFF
              IF(I.NE.NBUFF1)THEN
                IF(I.LE.NMAXBUFF/2)THEN
                  WRITE(CDUMMY,'(A2,I1,A1)') '#[',I,']'
                ELSE
                  WRITE(CDUMMY,'(A5,I1)') 'err #',I-NMAXBUFF/2
                END IF
                L=TRUELEN(CDUMMY)
                NB=10*(I-1)+NORIGEN+3
                CALL BUTTON(NB,CDUMMY(1:L),0)
                CALL BUTTON(NB,CDUMMY(1:L),3)
              END IF
            END DO
          END IF
C..............................................................................
C la operacion '<' puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'<')THEN
              NAXIS(1,NBUFF0)=NAXIS(1,NBUFF1)
              NAXIS(2,NBUFF0)=NAXIS(2,NBUFF1)
              NFRAMES(NBUFF0)=NFRAMES(NBUFF1)
              LWAVECAL(NBUFF0)=LWAVECAL(NBUFF1)
              CRPIX1(NBUFF0)=CRPIX1(NBUFF1)
              CRVAL1(NBUFF0)=CRVAL1(NBUFF1)
              CDELT1(NBUFF0)=CDELT1(NBUFF1)
              DO I=1,NAXIS(2,NBUFF0)
                DO J=1,NAXIS(1,NBUFF0)
                  IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF1)
                END DO
              END DO
              IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +           (NBUFF1.LE.NMAXBUFF/2))THEN
                NAXIS(1,NBUFF0+NMAXBUFF/2)=NAXIS(1,NBUFF1+NMAXBUFF/2)
                NAXIS(2,NBUFF0+NMAXBUFF/2)=NAXIS(2,NBUFF1+NMAXBUFF/2)
                NFRAMES(NBUFF0+NMAXBUFF/2)=NFRAMES(NBUFF1+NMAXBUFF/2)
                DO I=1,NAXIS(2,NBUFF0+NMAXBUFF/2)
                  DO J=1,NAXIS(1,NBUFF0+NMAXBUFF/2)
                    IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(J,I,NBUFF1+NMAXBUFF/2)
                  END DO
                END DO
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C la operacion '4' puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'4')THEN
              NAXIS(1,NBUFF0)=NAXIS(1,NBUFF1)
              NAXIS(2,NBUFF0)=NAXIS(2,NBUFF1)
              NFRAMES(NBUFF0)=NFRAMES(NBUFF1)
              NMED1=NAXIS(1,NBUFF0)/2
              NMED2=NAXIS(2,NBUFF0)/2
              DO I=1,NMED2
                DO J=1,NMED1
                  IMAGEN(J,I,NBUFF0)=IMAGEN(2*J-1,2*I-1,NBUFF1)
                  IMAGEN(J+NMED1,I,NBUFF0)=IMAGEN(2*J,2*I-1,NBUFF1)
                  IMAGEN(J,I+NMED2,NBUFF0)=IMAGEN(2*J-1,2*I,NBUFF1)
                  IMAGEN(J+NMED1,I+NMED2,NBUFF0)=IMAGEN(2*J,2*I,NBUFF1)
                END DO
              END DO
              IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +           (NBUFF1.LE.NMAXBUFF/2))THEN
                NAXIS(1,NBUFF0+NMAXBUFF/2)=NAXIS(1,NBUFF1+NMAXBUFF/2)
                NAXIS(2,NBUFF0+NMAXBUFF/2)=NAXIS(2,NBUFF1+NMAXBUFF/2)
                NFRAMES(NBUFF0+NMAXBUFF/2)=NFRAMES(NBUFF1+NMAXBUFF/2)
                NMED1=NAXIS(1,NBUFF0+NMAXBUFF/2)/2
                NMED2=NAXIS(2,NBUFF0+NMAXBUFF/2)/2
                DO I=1,NMED2
                  DO J=1,NMED1
                    IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(2*J-1,2*I-1,NBUFF1+NMAXBUFF/2)
                    IMAGEN(J+NMED1,I,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(2*J,2*I-1,NBUFF1+NMAXBUFF/2)
                    IMAGEN(J,I+NMED2,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(2*J-1,2*I,NBUFF1+NMAXBUFF/2)
                    IMAGEN(J+NMED1,I+NMED2,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(2*J,2*I,NBUFF1+NMAXBUFF/2)
                  END DO
                END DO
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C la operacion '1' puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'1')THEN
              NAXIS(1,NBUFF0)=NAXIS(1,NBUFF1)
              NAXIS(2,NBUFF0)=NAXIS(2,NBUFF1)
              NFRAMES(NBUFF0)=NFRAMES(NBUFF1)
              NMED1=NAXIS(1,NBUFF0)/2
              NMED2=NAXIS(2,NBUFF0)/2
              DO I=1,NMED2
                DO J=1,NMED1
                  IMAGEN(2*J-1,2*I-1,NBUFF0)=IMAGEN(J,I,NBUFF1)
                  IMAGEN(2*J,2*I-1,NBUFF0)=IMAGEN(J+NMED1,I,NBUFF1)
                  IMAGEN(2*J-1,2*I,NBUFF0)=IMAGEN(J,I+NMED2,NBUFF1)
                  IMAGEN(2*J,2*I,NBUFF0)=IMAGEN(J+NMED1,I+NMED2,NBUFF1)
                END DO
              END DO
              IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +           (NBUFF1.LE.NMAXBUFF/2))THEN
                NAXIS(1,NBUFF0+NMAXBUFF/2)=NAXIS(1,NBUFF1+NMAXBUFF/2)
                NAXIS(2,NBUFF0+NMAXBUFF/2)=NAXIS(2,NBUFF1+NMAXBUFF/2)
                NFRAMES(NBUFF0+NMAXBUFF/2)=NFRAMES(NBUFF1+NMAXBUFF/2)
                NMED1=NAXIS(1,NBUFF0+NMAXBUFF/2)/2
                NMED2=NAXIS(2,NBUFF0+NMAXBUFF/2)/2
                DO I=1,NMED2
                  DO J=1,NMED1
                    IMAGEN(2*J-1,2*I-1,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(J,I,NBUFF1+NMAXBUFF/2)
                    IMAGEN(2*J,2*I-1,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(J+NMED1,I,NBUFF1+NMAXBUFF/2)
                    IMAGEN(2*J-1,2*I,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(J,I+NMED2,NBUFF1+NMAXBUFF/2)
                    IMAGEN(2*J,2*I,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(J+NMED1,I+NMED2,NBUFF1+NMAXBUFF/2)
                  END DO
                END DO
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C la operacion '9' puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'9')THEN
              NAXIS(1,NBUFF0)=3*NAXIS(1,NBUFF1)
              NAXIS(2,NBUFF0)=3*NAXIS(2,NBUFF1)
              NFRAMES(NBUFF0)=9
              DO I=1,NAXIS(2,NBUFF1)
                DO J=1,NAXIS(1,NBUFF1)
                  DO K=1,9
                    IMAGEN(J+DJ(K),I+DI(K),NBUFF0)=IMAGEN(J,I,NBUFF1)
                  END DO
                END DO
              END DO
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C la operacion 'c' puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'c')THEN
              NAXIS(1,NBUFF0)=NAXIS(1,NBUFF1)/3
              NAXIS(2,NBUFF0)=NAXIS(2,NBUFF1)/3
              NFRAMES(NBUFF0)=1
              DO I=1,NAXIS(2,NBUFF1)
                DO J=1,NAXIS(1,NBUFF1)
                  IMAGEN(J,I,NBUFF0)=0.
                  DO K=1,9
                    IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)+
     +               IMAGEN(J+DJ(K),I+DI(K),NBUFF1)
                  END DO
                  IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)/9.0
                END DO
              END DO
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C la operacion 'q' puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'q')THEN
              NAXIS(1,NBUFF0)=NAXIS(1,NBUFF1)
              NAXIS(2,NBUFF0)=NAXIS(2,NBUFF1)
              NFRAMES(NBUFF0)=NFRAMES(NBUFF1)
              NMED1=NAXIS(1,NBUFF0)/2
              NMED2=NAXIS(2,NBUFF0)/2
              NQUAD0=READILIM('Quadrant no. in output buffer','@',1,4)
              IF(NQUAD0.EQ.1)THEN
                NOFFX0=512
                NOFFY0=512
              ELSEIF(NQUAD0.EQ.2)THEN
                NOFFX0=0
                NOFFY0=512
              ELSEIF(NQUAD0.EQ.3)THEN
                NOFFX0=0
                NOFFY0=0
              ELSE
                NOFFX0=512
                NOFFY0=0
              END IF
              NQUAD1=READILIM('Quadrant no. in  input buffer','@',1,4)
              IF(NQUAD1.EQ.1)THEN
                NOFFX1=512
                NOFFY1=512
              ELSEIF(NQUAD1.EQ.2)THEN
                NOFFX1=0
                NOFFY1=512
              ELSEIF(NQUAD1.EQ.3)THEN
                NOFFX1=0
                NOFFY1=0
              ELSE
                NOFFX1=512
                NOFFY1=0
              END IF
              DO I=1,NMED2
                DO J=1,NMED1
                  IMAGEN(J+NOFFX0,I+NOFFY0,NBUFF0)=
     +             IMAGEN(J+NOFFX1,I+NOFFY1,NBUFF1)
                END DO
              END DO
              IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +           (NBUFF1.LE.NMAXBUFF/2))THEN
                NAXIS(1,NBUFF0+NMAXBUFF/2)=NAXIS(1,NBUFF1+NMAXBUFF/2)
                NAXIS(2,NBUFF0+NMAXBUFF/2)=NAXIS(2,NBUFF1+NMAXBUFF/2)
                NFRAMES(NBUFF0+NMAXBUFF/2)=NFRAMES(NBUFF1+NMAXBUFF/2)
                NMED1=NAXIS(1,NBUFF0+NMAXBUFF/2)/2
                NMED2=NAXIS(2,NBUFF0+NMAXBUFF/2)/2
                DO I=1,NMED2
                  DO J=1,NMED1
                    IMAGEN(J+NOFFX0,I+NOFFY0,NBUFF0+NMAXBUFF/2)=
     +               IMAGEN(J+NOFFX1,I+NOFFY1,NBUFF1+NMAXBUFF/2)
                  END DO
                END DO
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C la operacion entre dos buffers (sin cortes en X o en Y) tambien puede hacerse
          IF(.NOT.LCANCEL)THEN
            IF(CUTIL.EQ.'f')THEN
              IF(COPER.EQ.'+')THEN
                IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +             (NBUFF1.LE.NMAXBUFF/2))THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)+
     +                 IMAGEN(J,I,NBUFF1+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF1+NMAXBUFF/2))
                    END DO
                  END DO
                END IF
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF0)=
     +               IMAGEN(J,I,NBUFF0)+IMAGEN(J,I,NBUFF1)
                  END DO
                END DO
              ELSEIF(COPER.EQ.'-')THEN
                IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +             (NBUFF1.LE.NMAXBUFF/2))THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)+
     +                 IMAGEN(J,I,NBUFF1+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF1+NMAXBUFF/2))
                    END DO
                  END DO
                END IF
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF0)=
     +               IMAGEN(J,I,NBUFF0)-IMAGEN(J,I,NBUFF1)
                  END DO
                END DO
              ELSEIF(COPER.EQ.'*')THEN
                IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +             (NBUFF1.LE.NMAXBUFF/2))THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                 IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                 IMAGEN(J,I,NBUFF1+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF1+NMAXBUFF/2)+
     +                 IMAGEN(J,I,NBUFF1)*IMAGEN(J,I,NBUFF1)*
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                 IMAGEN(J,I,NBUFF0+NMAXBUFF/2))
                    END DO
                  END DO
                END IF
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF0)=
     +               IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF1)
                  END DO
                END DO
              ELSEIF(COPER.EQ.'/')THEN
                LDIVZERO=.FALSE.
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IF(IMAGEN(J,I,NBUFF1).EQ.0.0) LDIVZERO=.TRUE.
                  END DO
                END DO
                IF(LDIVZERO)THEN
                  WRITE(*,101) '***ERROR***'
                  WRITE(*,101) '=> Division by zero!'
                  CZERO(1:1)=READC('e[x]it or [c]ontinue','c','xc')
                  IF(CZERO.EQ.'c')THEN
                    FZERO=READF('Pixel value for infinity','0.0')
                    IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +                 (NBUFF1.LE.NMAXBUFF/2))THEN
                      FZERO_ERR=
     +                 READF('Pixel value for error in infinity','0.0')
                      DO I=NY1,NY2
                        DO J=NX1,NX2
                          IF(IMAGEN(J,I,NBUFF1).EQ.0.0)THEN
                            IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=FZERO_ERR
                          ELSE
                            IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                       IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                       IMAGEN(J,I,NBUFF1+NMAXBUFF/2)*
     +                       IMAGEN(J,I,NBUFF1+NMAXBUFF/2)+
     +                       IMAGEN(J,I,NBUFF1)*IMAGEN(J,I,NBUFF1)*
     +                       IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                       IMAGEN(J,I,NBUFF0+NMAXBUFF/2))/
     +                       (IMAGEN(J,I,NBUFF1)*IMAGEN(J,I,NBUFF1))
                          END IF
                        END DO
                      END DO
                    END IF
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IF(IMAGEN(J,I,NBUFF1).EQ.0.0)THEN
                          IMAGEN(J,I,NBUFF0)=FZERO
                        ELSE
                          IMAGEN(J,I,NBUFF0)=
     +                     IMAGEN(J,I,NBUFF0)/IMAGEN(J,I,NBUFF1)
                        END IF
                      END DO
                    END DO
                  END IF
                ELSE
                  IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +               (NBUFF1.LE.NMAXBUFF/2))THEN
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                   IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                   IMAGEN(J,I,NBUFF1+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF1+NMAXBUFF/2)+
     +                   IMAGEN(J,I,NBUFF1)*IMAGEN(J,I,NBUFF1)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2))/
     +                   (IMAGEN(J,I,NBUFF1)*IMAGEN(J,I,NBUFF1))
                      END DO
                    END DO
                  END IF
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0)=
     +                 IMAGEN(J,I,NBUFF0)/IMAGEN(J,I,NBUFF1)
                    END DO
                  END DO
                END IF
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C la operacion 'f' (filter) puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'f')THEN
              WRITE(*,101) ' '
              WRITE(*,101) '(0) Compute local r.m.s.'
              WRITE(*,101) '(1) 2-D Mean filter'
              WRITE(*,101) '(2) 2-D Mean filter excluding points'
              WRITE(*,101) '(3) 2-D Median filter'
              WRITE(*,101) '(4) 1-D Savitzky-Golay filter'
              WRITE(*,100) '(5) 1-D Savitzky-Golay filter'
              WRITE(*,101) ' (excluding points)'
              WRITE(*,101) '(6) 1-D: 0.5,1,...,1,...,1,0.5 filter'
              WRITE(*,101) '(7) 1-D: 1,0,1,0,...,1,...,0,1,0,1 filter'
              WRITE(*,101) '(8) Set to zero below threshold'
              WRITE(*,101) '(9) Set to constant for data in a range'
              WRITE(*,101) '(k) Apply predefined kernel'
              WRITE(*,101) '(m) Mathematical transformation'
              WRITE(*,101) '(x) Reverse image in the x direction'
              WRITE(*,101) '(y) Reverse image in the y direction'
              WRITE(*,*)
              WRITE(*,100) '* Note: so far these operations are not'
              WRITE(*,101) ' translated to the error frames,'
              WRITE(*,100) '  except for option (2), in which the '
              WRITE(*,101) 'r.m.s. is stored in the associated '
              WRITE(*,101) '  error buffer'
              WRITE(*,*)
              CFILT(1:1)=READC('Option',CFILT,'0123456789kmxy')
c
              NAXIS(1,NBUFF1)=NAXIS(1,NBUFF0)
              NAXIS(2,NBUFF1)=NAXIS(2,NBUFF0)
              NFRAMES(NBUFF1)=NFRAMES(NBUFF0)
              DO I=1,NAXIS(2,NBUFF1)
                DO J=1,NAXIS(1,NBUFF1)
                  IMAGEN(J,I,NBUFF1)=0.
                END DO
              END DO
              IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +           (NBUFF1.LE.NMAXBUFF/2))THEN
                NAXIS(1,NBUFF1+NMAXBUFF/2)=NAXIS(1,NBUFF0+NMAXBUFF/2)
                NAXIS(2,NBUFF1+NMAXBUFF/2)=NAXIS(2,NBUFF0+NMAXBUFF/2)
                NFRAMES(NBUFF1+NMAXBUFF/2)=NFRAMES(NBUFF0+NMAXBUFF/2)
                DO I=1,NAXIS(2,NBUFF1+NMAXBUFF/2)
                  DO J=1,NAXIS(1,NBUFF1+NMAXBUFF/2)
                    IMAGEN(J,I,NBUFF1+NMAXBUFF/2)=0.
                  END DO
                END DO
              END IF
c set options for selected filter
              IF(CFILT.EQ.'0')THEN
                WRITE(*,101) '(1) traditional r.m.s.'
                WRITE(*,101) '(2) traditional r.m.s. excluding points'
                WRITE(*,101) '(3) robust r.m.s.'
                WRITE(*,101) '(4) robust r.m.s. excluding points'
                CRMS(1:1)=READC('Option',CRMS,'1234')
                IF((CRMS.EQ.'2').OR.(CRMS.EQ.'4'))THEN
                  WRITE(CDUMMY,*) TSIGMA
                  TSIGMA=READF('Times sigma to exclude points',CDUMMY)
                END IF
              ELSEIF(CFILT.EQ.'1')THEN
              ELSEIF(CFILT.EQ.'2')THEN
                WRITE(CDUMMY,*) TSIGMA
                TSIGMA=READF('Times sigma to exclude points',CDUMMY)
              ELSEIF(CFILT.EQ.'3')THEN
              ELSEIF((CFILT.EQ.'4').OR.(CFILT.EQ.'5').OR.
     +         (CFILT.EQ.'6').OR.(CFILT.EQ.'7'))THEN
                CAXIS(1:1)=READC('Axis direction (x/y)',CAXIS,'xy')
                WRITE(CDUMMY,*) IWIDTH
                IF(CAXIS.EQ.'x')THEN
                  IWIDTH=READILIM('Window width (pixels; odd)',CDUMMY,
     +             3,NAXIS(1,NBUFF0))
                ELSE
                  IWIDTH=READILIM('Window width (pixels; odd)',CDUMMY,
     +             3,NAXIS(2,NBUFF0))
                END IF
                IF(MOD(IWIDTH,2).EQ.0)THEN
                  IF(IWIDTH.GE.4)THEN
                    IWIDTH=IWIDTH-1
                  ELSE
                    IWIDTH=IWIDTH+1
                  END IF
                  WRITE(*,101) '***WARNING***'
                  WRITE(*,100) 'Window width must be odd => New width: '
                  WRITE(*,*) IWIDTH
                END IF
                IF((CFILT.EQ.'4').OR.(CFILT.EQ.'5'))THEN
                  WRITE(CDUMMY,*) NDEGREE
                  NDEGREE=READILIM('Polynomial degree',CDUMMY,
     +             0,MIN0(19,IWIDTH-1))
                  IF(CFILT.EQ.'5')THEN
                    WRITE(CDUMMY,*) TSIGMA
                    TSIGMA=READF('Times sigma to exclude points',CDUMMY)
                  END IF
                  CFUN_OR_DER(1:1)=READC('Function (1), derivative (2)',
     +             CFUN_OR_DER,'12')
                END IF
              ELSEIF(CFILT.EQ.'8')THEN
                THRESHOLD=READF('Threshold','@')
              ELSEIF(CFILT.EQ.'9')THEN
                NEWCONSTANT=READF('New constant','@')
                THRESHOLD1=READF('Mininum value for range','@')
                THRESHOLD2=READF('Maximum value for range','@')
              END IF
c execute filter
              NEXTINFO=0
              IF(CFILT.EQ.'x')THEN
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF1)=IMAGEN(NX2+NX1-J,I,NBUFF0)
                  END DO
                END DO
              ELSEIF(CFILT.EQ.'y')THEN
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF1)=IMAGEN(J,NY2+NY1-I,NBUFF0)
                  END DO
                END DO
              ELSEIF(CFILT.EQ.'m')THEN
                WRITE(*,101)'* Only logarithmic transformation so far!'
                OFFSET_LOG=READF('Offset for logarithmic transformation'
     +           , '@')
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IMAGEN(J,I,NBUFF1)=
     +               ALOG10(IMAGEN(J,I,NBUFF0)+OFFSET_LOG)
                  END DO
                END DO
              ELSEIF(CFILT.EQ.'k')THEN
                WRITE(*,101)'* Only sharpen kernel so far!'
                CKERNEL=READC('kernel side (3, 5, 7 o 9)','3','3579')
                IF(CKERNEL.EQ.'3')THEN
                  IKERNEL=1
                ELSEIF(CKERNEL.EQ.'5')THEN
                  IKERNEL=2
                ELSEIF(CKERNEL.EQ.'7')THEN
                  IKERNEL=3
                ELSE
                  IKERNEL=4
                END IF
                NX1_=NX1+IKERNEL
                NX2_=NX2-IKERNEL
                NY1_=NY1+IKERNEL
                NY2_=NY2-IKERNEL
                DO I=NY1_,NY2_
                  DO J=NX1_,NX2_
                    FDUMMY=0.0
                    DO II=I-IKERNEL,I+IKERNEL
                      DO JJ=J-IKERNEL,J+IKERNEL
                        IF((II.EQ.I).AND.(JJ.EQ.J))THEN
                          FDUMMY=FDUMMY+FLOAT(2*IKERNEL+1)
     +                     *FLOAT(2*IKERNEL+1)*IMAGEN(JJ,II,NBUFF0)
                        ELSE
                          FDUMMY=FDUMMY-IMAGEN(JJ,II,NBUFF0)
                        END IF
                      END DO
                    END DO
                    IMAGEN(J,I,NBUFF1)=FDUMMY
                  END DO
                END DO
              ELSEIF(CFILT.EQ.'8')THEN
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IF(IMAGEN(J,I,NBUFF0).LT.THRESHOLD)THEN
                      IMAGEN(J,I,NBUFF1)=0.
                    ELSE
                      IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF0)
                    END IF
                  END DO
                END DO
              ELSEIF(CFILT.EQ.'9')THEN
                DO I=NY1,NY2
                  DO J=NX1,NX2
                    IF((IMAGEN(J,I,NBUFF0).GE.THRESHOLD1).AND.
     +                 (IMAGEN(J,I,NBUFF0).LE.THRESHOLD2))THEN
                      IMAGEN(J,I,NBUFF1)=NEWCONSTANT
                    ELSE
                      IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF0)
                    END IF
                  END DO
                END DO
              ELSEIF(CFILT.EQ.'7')THEN
                IF(CAXIS.EQ.'x')THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      NX1_=J-IWIDTH/2
                      NX2_=J+IWIDTH/2
                      IF(NX1_.LT.NX1) NX1_=NX1
                      IF(NX2_.GT.NX2) NX2_=NX2
                      IMAGEN(J,I,NBUFF1)=0.
                      SUMFACTOR=0.
                      DO JJ=NX1_,NX2_
                        IF(MOD(JJ-J,2).EQ.0)THEN
                          FACTOR=1.
                        ELSE
                          FACTOR=0.
                        END IF
                        IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF1)+
     +                   IMAGEN(JJ,I,NBUFF0)*FACTOR
                        SUMFACTOR=SUMFACTOR+FACTOR
                      END DO
                      IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF1)/SUMFACTOR
                    END DO
                    CALL SHOWPERC(NY1,NY2,1,I,NEXTINFO)
                  END DO
                ELSE
                  DO J=NX1,NX2
                    DO I=NY1,NY2
                      NY1_=I-IWIDTH/2
                      NY2_=I+IWIDTH/2
                      IF(NY1_.LT.NY1) NY1_=NY1
                      IF(NY2_.GT.NY2) NY2_=NY2
                      IMAGEN(J,I,NBUFF1)=0.
                      SUMFACTOR=0.
                      DO II=NY1_,NY2_
                        IF(MOD(II-I,2).EQ.0)THEN
                          FACTOR=1.
                        ELSE
                          FACTOR=0.
                        END IF
                        IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF1)+
     +                   IMAGEN(J,II,NBUFF0)*FACTOR
                        SUMFACTOR=SUMFACTOR+FACTOR
                      END DO
                      IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF1)/SUMFACTOR
                    END DO
                    CALL SHOWPERC(NX1,NX2,1,J,NEXTINFO)
                  END DO
                END IF
              ELSEIF(CFILT.EQ.'6')THEN
                IF(CAXIS.EQ.'x')THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      NX1_=J-IWIDTH/2
                      NX2_=J+IWIDTH/2
                      IF(NX1_.LT.NX1) NX1_=NX1
                      IF(NX2_.GT.NX2) NX2_=NX2
                      IMAGEN(J,I,NBUFF1)=0.
                      SUMFACTOR=0.
                      DO JJ=NX1_,NX2_
                        IF(JJ.EQ.J-IWIDTH/2)THEN
                          FACTOR=0.5
                        ELSEIF(JJ.EQ.J+IWIDTH/2)THEN
                          FACTOR=0.5
                        ELSE
                          FACTOR=1.
                        END IF
                        IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF1)+
     +                   IMAGEN(JJ,I,NBUFF0)*FACTOR
                        SUMFACTOR=SUMFACTOR+FACTOR
                      END DO
                      IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF1)/SUMFACTOR
                    END DO
                    CALL SHOWPERC(NY1,NY2,1,I,NEXTINFO)
                  END DO
                ELSE
                  DO J=NX1,NX2
                    DO I=NY1,NY2
                      NY1_=I-IWIDTH/2
                      NY2_=I+IWIDTH/2
                      IF(NY1_.LT.NY1) NY1_=NY1
                      IF(NY2_.GT.NY2) NY2_=NY2
                      IMAGEN(J,I,NBUFF1)=0.
                      SUMFACTOR=0.
                      DO II=NY1_,NY2_
                        IF(II.EQ.I-IWIDTH/2)THEN
                          FACTOR=0.5
                        ELSEIF(II.EQ.I+IWIDTH/2)THEN
                          FACTOR=0.5
                        ELSE
                          FACTOR=1.
                        END IF
                        IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF1)+
     +                   IMAGEN(J,II,NBUFF0)*FACTOR
                        SUMFACTOR=SUMFACTOR+FACTOR
                      END DO
                      IMAGEN(J,I,NBUFF1)=IMAGEN(J,I,NBUFF1)/SUMFACTOR
                    END DO
                    CALL SHOWPERC(NX1,NX2,1,J,NEXTINFO)
                  END DO
                END IF
              ELSEIF((CFILT.EQ.'4').OR.(CFILT.EQ.'5'))THEN
                ISYSTEM=SYSTEMFUNCTION('date')
                IF(CAXIS.EQ.'x')THEN
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      NX1_=J-IWIDTH/2
                      NX2_=J+IWIDTH/2
                      IF(NX1_.LT.NX1) NX1_=NX1
                      IF(NX2_.GT.NX2) NX2_=NX2
                      K=0
                      DO JJ=NX1_,NX2_
                        K=K+1
                        XF(K)=REAL(JJ)
                        YF(K)=IMAGEN(JJ,I,NBUFF0)
                      END DO
                      IF(K.GE.NDEGREE+1)THEN
                        NCOEF=NDEGREE+1
                      ELSE
                        NCOEF=K-1
                      END IF
                      IF(CFILT.EQ.'4')THEN
                        CALL POLFIT(XF,YF,YF,K,NCOEF,0,COEF,CHISQR)
                      ELSE
                        CALL POLFITSIG(K,XF,YF,TSIGMA,NCOEF-1,COEF,LFIT)
                      END IF
                      IF(CFUN_OR_DER.EQ.'1')THEN !function
                        IMAGEN(J,I,NBUFF1)=COEF(NCOEF)
                        IF(NCOEF.GT.1)THEN
                          DO K=NCOEF-1,1,-1
                            IMAGEN(J,I,NBUFF1)=
     +                       IMAGEN(J,I,NBUFF1)*REAL(J)+COEF(K)
                          END DO
                        END IF
                      ELSE                       !derivative
                        IMAGEN(J,I,NBUFF1)=0.0
                        IF(NCOEF.GT.1)THEN
                          DO K=2,NCOEF
                            COEFDER(K-1)=(K-1)*COEF(K)
                          END DO
                          NCOEFDER=NCOEF-1
                          IMAGEN(J,I,NBUFF1)=COEFDER(NCOEFDER)
                          IF(NCOEFDER.GT.1)THEN
                            DO K=NCOEFDER-1,1,-1
                              IMAGEN(J,I,NBUFF1)=
     +                         IMAGEN(J,I,NBUFF1)*REAL(J)+COEFDER(K)
                            END DO
                          END IF
                        END IF
                      END IF
                    END DO
                    CALL SHOWPERC(NY1,NY2,1,I,NEXTINFO)
                  END DO
                ELSE
                  DO J=NX1,NX2
                    DO I=NY1,NY2
                      NY1_=I-IWIDTH/2
                      NY2_=I+IWIDTH/2
                      IF(NY1_.LT.NY1) NY1_=NY1
                      IF(NY2_.GT.NY2) NY2_=NY2
                      K=0
                      DO II=NY1_,NY2_
                        K=K+1
                        XF(K)=REAL(II)
                        YF(K)=IMAGEN(J,II,NBUFF0)
                      END DO
                      IF(K.GE.NDEGREE+1)THEN
                        NCOEF=NDEGREE+1
                      ELSE
                        NCOEF=K-1
                      END IF
                      IF(CFILT.EQ.'4')THEN
                        CALL POLFIT(XF,YF,YF,K,NCOEF,0,COEF,CHISQR)
                      ELSE
                        CALL POLFITSIG(K,XF,YF,TSIGMA,NCOEF-1,COEF,LFIT)
                      END IF
                      IF(CFUN_OR_DER.EQ.'1')THEN !function
                        IMAGEN(J,I,NBUFF1)=COEF(NCOEF)
                        IF(NCOEF.GT.1)THEN
                          DO K=NCOEF-1,1,-1
                            IMAGEN(J,I,NBUFF1)=
     +                       IMAGEN(J,I,NBUFF1)*REAL(I)+COEF(K)
                          END DO
                        END IF
                      ELSE                       !derivative
                        IMAGEN(J,I,NBUFF1)=0.0
                        IF(NCOEF.GT.1)THEN
                          DO K=2,NCOEF
                            COEFDER(K-1)=(K-1)*COEF(K)
                          END DO
                          NCOEFDER=NCOEF-1
                          IMAGEN(J,I,NBUFF1)=COEFDER(NCOEFDER)
                          IF(NCOEFDER.GT.1)THEN
                            DO K=NCOEFDER-1,1,-1
                              IMAGEN(J,I,NBUFF1)=
     +                         IMAGEN(J,I,NBUFF1)*REAL(I)+COEFDER(K)
                            END DO
                          END IF
                        END IF
                      END IF
                    END DO
                    CALL SHOWPERC(NX1,NX2,1,J,NEXTINFO)
                  END DO
                END IF
                ISYSTEM=SYSTEMFUNCTION('date')
              ELSE ! CFILT = 0, 1, 2 or 3
                WRITE(CDUMMY,*) NBOX_X
                NBOX_X=READILIM('Box width in X (must be odd)',
     +           CDUMMY,1,NAXIS(1,NBUFF0))
                IF(MOD(NBOX_X,2).EQ.0)THEN
                  IF(NBOX_X.GE.4)THEN
                    NBOX_X=NBOX_X-1
                  ELSE
                    NBOX_X=NBOX_X+1
                  END IF
                  WRITE(*,101) '***WARNING***'
                  WRITE(*,100) 'Window width must be odd => New width: '
                  WRITE(*,*) NBOX_X
                END IF
                WRITE(CDUMMY,*) NBOX_Y
                NBOX_Y=READILIM('Box width in Y (must be odd)',
     +           CDUMMY,1,NAXIS(2,NBUFF0))
                IF(MOD(NBOX_Y,2).EQ.0)THEN
                  IF(NBOX_Y.GE.4)THEN
                    NBOX_Y=NBOX_Y-1
                  ELSE
                    NBOX_Y=NBOX_Y+1
                  END IF
                  WRITE(*,101) '***WARNING***'
                  WRITE(*,100) 'Window width must be odd => New width: '
                  WRITE(*,*) NBOX_Y
                END IF
                ISYSTEM=SYSTEMFUNCTION('date')
                DO I=NY1,NY2
                  NY1_=I-NBOX_Y/2
                  NY2_=I+NBOX_Y/2
                  IF(NY1_.LT.NY1) NY1_=NY1
                  IF(NY2_.GT.NY2) NY2_=NY2
                  DO J=NX1,NX2
                    NX1_=J-NBOX_X/2
                    NX2_=J+NBOX_X/2
                    IF(NX1_.LT.NX1) NX1_=NX1
                    IF(NX2_.GT.NX2) NX2_=NX2
                    K=0
                    DO II=NY1_,NY2_
                      DO JJ=NX1_,NX2_
                        K=K+1
                        PIXEL(K)=IMAGEN(JJ,II,NBUFF0)
                      END DO
                    END DO
                    IF(CFILT.EQ.'0')THEN
                      IF(CRMS.EQ.'1')THEN
                        FDUMMY=FMEAN0(K,PIXEL,FSIGMA)
                      ELSEIF(CRMS.EQ.'2')THEN
                        FDUMMY=FMEAN2(K,PIXEL,TSIGMA,FSIGMA)
                      ELSEIF(CRMS.EQ.'3')THEN
                        FSIGMA=FROBUSTRMS0(K,PIXEL)
                      ELSEIF(CRMS.EQ.'4')THEN
                        FSIGMA=FROBUSTRMS1(K,PIXEL,TSIGMA)
                      ELSE
                        WRITE(*,101) '***FATAL ERROR***'
                        WRITE(*,101) '=> invalid CRMS value: '//CRMS
                        STOP
                      END IF
                      IMAGEN(J,I,NBUFF1)=FSIGMA
                    ELSEIF(CFILT.EQ.'1')THEN
                      IMAGEN(J,I,NBUFF1)=FMEAN1(K,PIXEL)
                    ELSEIF(CFILT.EQ.'2')THEN
                      IMAGEN(J,I,NBUFF1)=FMEAN2(K,PIXEL,TSIGMA,FSIGMA)
                      IMAGEN(J,I,NBUFF1+NMAXBUFF/2)=FSIGMA
                    ELSEIF(CFILT.EQ.'3')THEN
                      IMAGEN(J,I,NBUFF1)=FMEDIAN1(K,PIXEL)
                    ELSE
                      WRITE(*,101) '***FATAL ERROR***'
                      WRITE(*,101) '=> invalid CFILT value: '//CFILT
                      STOP
                    END IF
                  END DO
                  CALL SHOWPERC(NY1,NY2,1,I,NEXTINFO)
                END DO
                ISYSTEM=SYSTEMFUNCTION('date')
              END IF
c
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C la operacion 's' (offset xy) puede hacerse ya
          IF(.NOT.LCANCEL)THEN
            IF(COPER.EQ.'s')THEN
              XOFFSET=READF('Offset for X shift (+right,-left)',
     +         '0.0')
              YOFFSET=READF('Offset for Y shift ...(+up,-down)',
     +         '0.0')
              J0=INT(XOFFSET)
              FJ0=XOFFSET-REAL(J0)
              I0=INT(YOFFSET)
              FI0=YOFFSET-REAL(I0)
c inicializamos imagen de destino a cero
              NAXIS(1,NBUFF1)=NAXIS(1,NBUFF0)
              NAXIS(2,NBUFF1)=NAXIS(2,NBUFF0)
              NFRAMES(NBUFF1)=NFRAMES(NBUFF0)
              DO I=1,NAXIS(2,NBUFF1)
                DO J=1,NAXIS(1,NBUFF1)
                  IMAGEN(J,I,NBUFF1)=0.
                END DO
              END DO
c si tenemos errores, inicializamos imagen de errores de destino a cero
              IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +           (NBUFF1.LE.NMAXBUFF/2))THEN
                NAXIS(1,NBUFF1+NMAXBUFF/2)=NAXIS(1,NBUFF0+NMAXBUFF/2)
                NAXIS(2,NBUFF1+NMAXBUFF/2)=NAXIS(2,NBUFF0+NMAXBUFF/2)
                NFRAMES(NBUFF1+NMAXBUFF/2)=NFRAMES(NBUFF0+NMAXBUFF/2)
                DO I=1,NAXIS(2,NBUFF1+NMAXBUFF/2)
                  DO J=1,NAXIS(1,NBUFF1+NMAXBUFF/2)
                    IMAGEN(J,I,NBUFF1+NMAXBUFF/2)=0.
                  END DO
                END DO
              END IF
c aplicamos el offset a la imagen de datos (I,J son el pixel en la imagen
C inicial, mientras que I1,J y I2,J son los dos posibles pixels a los que
C movemos la sen~al en la nueva imagen)
              DO I=1,NAXIS(2,NBUFF0)
                I1=I+I0
                IF(FI0.EQ.0.0)THEN
                  I2=0
                ELSEIF(FI0.GT.0.0)THEN
                  I2=I1+1
                ELSEIF(FI0.LT.0.0)THEN
                  I2=I1-1
                END IF
                DO J=1,NAXIS(1,NBUFF0)
                  J1=J+J0
                  IF(FJ0.EQ.0.0)THEN
                    J2=0
                  ELSEIF(FJ0.GT.0.0)THEN
                    J2=J1+1
                  ELSEIF(FJ0.LT.0.0)THEN
                    J2=J1-1
                  END IF
                  IF((J1.GE.1).AND.(J1.LE.NAXIS(1,NBUFF1)))THEN
                    IF((I1.GE.1).AND.(I1.LE.NAXIS(2,NBUFF1)))THEN
                      IMAGEN(J1,I1,NBUFF1)=IMAGEN(J1,I1,NBUFF1)+
     +                 IMAGEN(J,I,NBUFF0)*(1.-ABS(FI0))*(1.-ABS(FJ0))
                    END IF
                    IF((I2.GE.1).AND.(I2.LE.NAXIS(2,NBUFF1)))THEN
                      IMAGEN(J1,I2,NBUFF1)=IMAGEN(J1,I2,NBUFF1)+
     +                 IMAGEN(J,I,NBUFF0)*ABS(FI0)*(1.-ABS(FJ0))
                    END IF
                  END IF
                  IF((J2.GE.1).AND.(J2.LE.NAXIS(1,NBUFF1)))THEN
                    IF((I1.GE.1).AND.(I1.LE.NAXIS(2,NBUFF1)))THEN
                      IMAGEN(J2,I1,NBUFF1)=IMAGEN(J2,I1,NBUFF1)+
     +                 IMAGEN(J,I,NBUFF0)*(1.-ABS(FI0))*ABS(FJ0)
                    END IF
                    IF((I2.GE.1).AND.(I2.LE.NAXIS(2,NBUFF1)))THEN
                      IMAGEN(J2,I2,NBUFF1)=IMAGEN(J2,I2,NBUFF1)+
     +                 IMAGEN(J,I,NBUFF0)*ABS(FI0)*ABS(FJ0)
                    END IF
                  END IF
                END DO
              END DO
c si tenemos errores, aplicamos el offset a la imagen de errores
              IF((NBUFF0.LE.NMAXBUFF/2).AND.
     +           (NBUFF1.LE.NMAXBUFF/2))THEN
                DO I=1,NAXIS(2,NBUFF0+NMAXBUFF/2)
                  I1=I+I0
                  IF(FI0.EQ.0.0)THEN
                    I2=0
                  ELSEIF(FI0.GT.0.0)THEN
                    I2=I1+1
                  ELSEIF(FI0.LT.0.0)THEN
                    I2=I1-1
                  END IF
                  DO J=1,NAXIS(1,NBUFF0+NMAXBUFF/2)
                    J1=J+J0
                    IF(FJ0.EQ.0.0)THEN
                      J2=0
                    ELSEIF(FJ0.GT.0.0)THEN
                      J2=J1+1
                    ELSEIF(FJ0.LT.0.0)THEN
                      J2=J1-1
                    END IF
                    IF((J1.GE.1).AND.(J1.LE.NAXIS(1,NBUFF1)))THEN
                      IF((I1.GE.1).AND.
     +                 (I1.LE.NAXIS(2,NBUFF1+NMAXBUFF/2)))THEN
                        IMAGEN(J1,I1,NBUFF1+NMAXBUFF/2)=
     +                   IMAGEN(J1,I1,NBUFF1+NMAXBUFF/2)+
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   (1.-ABS(FI0))*(1.-ABS(FJ0))
                      END IF
                      IF((I2.GE.1).AND.
     +                 (I2.LE.NAXIS(2,NBUFF1+NMAXBUFF/2)))THEN
                        IMAGEN(J1,I2,NBUFF1+NMAXBUFF/2)=
     +                   IMAGEN(J1,I2,NBUFF1+NMAXBUFF/2)+
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   ABS(FI0)*(1.-ABS(FJ0))
                      END IF
                    END IF
                    IF((J2.GE.1).AND.(J2.LE.NAXIS(1,NBUFF1)))THEN
                      IF((I1.GE.1).AND.
     +                 (I1.LE.NAXIS(2,NBUFF1+NMAXBUFF/2)))THEN
                        IMAGEN(J2,I1,NBUFF1+NMAXBUFF/2)=
     +                   IMAGEN(J2,I1,NBUFF1+NMAXBUFF/2)+
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   (1.-ABS(FI0))*ABS(FJ0)
                      END IF
                      IF((I2.GE.1).AND.
     +                 (I2.LE.NAXIS(2,NBUFF1+NMAXBUFF/2)))THEN
                        IMAGEN(J2,I2,NBUFF1+NMAXBUFF/2)=
     +                   IMAGEN(J2,I2,NBUFF1+NMAXBUFF/2)+
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   ABS(FI0)*ABS(FJ0)
                      END IF
                    END IF
                  END DO
                END DO
                DO I=1,NAXIS(2,NBUFF1+NMAXBUFF/2)
                  DO J=1,NAXIS(1,NBUFF1+NMAXBUFF/2)
                    IMAGEN(J,I,NBUFF1+NMAXBUFF/2)=
     +               SQRT(IMAGEN(J,I,NBUFF1+NMAXBUFF/2))
                  END DO
                END DO
              END IF
c ya hemos terminado
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C operacion con los cortes
          IF(.NOT.LCANCEL)THEN
            WRITE(*,101) '(1*) mean cut'
            WRITE(*,101) '(2*) mean cut (+/- 3 sigma)'
            WRITE(*,101) '(3) median cut'
            WRITE(*,101) '(0) exit'
            ICUT=READILIM('Option (note: * with errors)','0',0,3)
            LCANCEL=(ICUT.EQ.0)
          END IF
C..............................................................................
C calculamos corte en X
          IF(.NOT.LCANCEL)THEN
            IF(CUTIL.EQ.'x')THEN
              DO I=1,NAXIS(2,NBUFF1)
                IFSCAN(I)=.FALSE.
              END DO
              LASK=.TRUE.
              DO WHILE(LASK)
                WRITE(CDUMMY,*) NY1_PLOT
                NY1_=READILIM('Y min (pixel, 0=EXIT)',CDUMMY,
     +           0,NAXIS(2,NBUFF1))
                LASK=(NY1_.NE.0)
                IF(LASK)THEN
                  WRITE(CDUMMY,*) NY2_PLOT
                  NY2_=READILIM('Y max (pixel)',CDUMMY,
     +             NY1_,NAXIS(2,NBUFF1))
                  LASK=(NY2_.NE.0)
                  IF(LASK)THEN
                    DO I=NY1_,NY2_
                      IFSCAN(I)=.TRUE.
                    END DO
                  END IF
                END IF
              END DO
              K=0
              DO I=1,NAXIS(2,NBUFF1)
                IF(IFSCAN(I)) K=K+1
              END DO
              LCANCEL=(K.EQ.0)
              IF(.NOT.LCANCEL)THEN
                DO J=1,NAXIS(1,NBUFF1)
                  K=0
                  DO I=1,NAXIS(2,NBUFF1)
                    IF(IFSCAN(I))THEN
                      K=K+1
                      YCUT(K)=IMAGEN(J,I,NBUFF1)
                      IF(NBUFF1.LE.NMAXBUFF/2)THEN
                        YCUT_ERR(K)=IMAGEN(J,I,NBUFF1+NMAXBUFF/2)
                      ELSE
                        YCUT_ERR(K)=0.
                      END IF
                    END IF
                  END DO
                  IF(ICUT.EQ.1)THEN
                    XCUT(J)=FMEAN0E(K,YCUT,YCUT_ERR,FSIGMA,XCUT_ERR(J))
                  ELSEIF(ICUT.EQ.2)THEN
                    XCUT(J)=
     +               FMEAN2E(K,YCUT,YCUT_ERR,3.0,FSIGMA,XCUT_ERR(J))
                  ELSEIF(ICUT.EQ.3)THEN
                    XCUT(J)=FMEDIAN1(K,YCUT)
                    XCUT_ERR(J)=0.0
                  END IF
                END DO
                IF(COPER.EQ.'+')THEN
                  IF(NBUFF0.LE.NMAXBUFF/2)THEN
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)+
     +                   XCUT_ERR(J)*XCUT_ERR(J))
                      END DO
                    END DO
                  END IF
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)+XCUT(J)
                    END DO
                  END DO
                ELSEIF(COPER.EQ.'-')THEN
                  IF(NBUFF0.LE.NMAXBUFF/2)THEN
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)+
     +                   XCUT_ERR(J)*XCUT_ERR(J))
                      END DO
                    END DO
                  END IF
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)-XCUT(J)
                    END DO
                  END DO
                ELSEIF(COPER.EQ.'*')THEN
                  IF(NBUFF0.LE.NMAXBUFF/2)THEN
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   XCUT(J)*XCUT(J)+
     +                   IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                   XCUT_ERR(J)*XCUT_ERR(J))
                      END DO
                    END DO
                  END IF
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)*XCUT(J)
                    END DO
                  END DO
                ELSEIF(COPER.EQ.'/')THEN
                  LDIVZERO=.FALSE.
                  DO J=NX1,NX2
                    IF(XCUT(J).EQ.0.0) LDIVZERO=.TRUE.
                  END DO
                  IF(LDIVZERO)THEN
                    WRITE(*,101) '***ERROR***'
                    WRITE(*,101) '=> Division by zero!'
                    CZERO(1:1)=READC('e[x]it or [c]ontinue','c','xc')
                    IF(CZERO.EQ.'c')THEN
                      FZERO=READF('Pixel value for infinity','0.0')
                      IF(NBUFF0.LE.NMAXBUFF/2)THEN
                        FZERO_ERR=READF(
     +                   'Pixel value for error at infinity','0.0')
                        IF(XCUT(J).EQ.0.0)THEN
                          DO I=NY1,NY2
                            DO J=NX1,NX2
                              IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=FZERO_ERR
                            END DO
                          END DO
                        ELSE
                          DO I=NY1,NY2
                            DO J=NX1,NX2
                              IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                         IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                         IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                         XCUT(J)*XCUT(J)+
     +                         IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                         XCUT_ERR(J)*XCUT_ERR(J))/
     +                         (XCUT(J)*XCUT(J))
                            END DO
                          END DO
                        END IF
                      END IF
                      DO J=NX1,NX2
                        IF(XCUT(J).EQ.0.0)THEN
                          DO I=NY1,NY2
                            IMAGEN(J,I,NBUFF0)=FZERO
                          END DO
                        ELSE
                          DO I=NY1,NY2
                            IMAGEN(J,I,NBUFF0)=
     +                       IMAGEN(J,I,NBUFF0)/XCUT(J)
                          END DO
                        END IF
                      END DO
                    END IF
                  ELSE
                    IF(NBUFF0.LE.NMAXBUFF/2)THEN
                      DO I=NY1,NY2
                        DO J=NX1,NX2
                          IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                     IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                     IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                     XCUT(J)*XCUT(J)+
     +                     IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                     XCUT_ERR(J)*XCUT_ERR(J))/
     +                     (XCUT(J)*XCUT(J))
                        END DO
                      END DO
                    END IF
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0)=
     +                   IMAGEN(J,I,NBUFF0)/XCUT(J)
                      END DO
                    END DO
                  END IF
                END IF
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
C calculamos corte en Y
          IF(.NOT.LCANCEL)THEN
            IF(CUTIL.EQ.'y')THEN
              DO J=1,NAXIS(1,NBUFF1)
                IFCHAN(J)=.FALSE.
              END DO
              LASK=.TRUE.
              DO WHILE(LASK)
                WRITE(CDUMMY,*) NX1_PLOT
                NX1_=READILIM('X min (pixel, 0=EXIT)',CDUMMY,
     +           0,NAXIS(1,NBUFF1))
                LASK=(NX1_.NE.0)
                IF(LASK)THEN
                  WRITE(CDUMMY,*) NX2_PLOT
                  NX2_=READILIM('X max (pixel)',CDUMMY,
     +             NX1_,NAXIS(1,NBUFF1))
                  LASK=(NX2_.NE.0)
                  IF(LASK)THEN
                    DO J=NX1_,NX2_
                      IFCHAN(J)=.TRUE.
                    END DO
                  END IF
                END IF
              END DO
              K=0
              DO J=1,NAXIS(1,NBUFF1)
                IF(IFCHAN(J)) K=K+1
              END DO
              LCANCEL=(K.EQ.0)
              IF(.NOT.LCANCEL)THEN
                DO I=1,NAXIS(2,NBUFF1)
                  K=0
                  DO J=1,NAXIS(1,NBUFF1)
                    IF(IFCHAN(J))THEN
                      K=K+1
                      XCUT(K)=IMAGEN(J,I,NBUFF1)
                      IF(NBUFF1.LE.NMAXBUFF/2)THEN
                        XCUT_ERR(K)=IMAGEN(J,I,NBUFF1+NMAXBUFF/2)
                      ELSE
                        XCUT_ERR(K)=0.
                      END IF
                    END IF
                  END DO
                  IF(ICUT.EQ.1)THEN
                    YCUT(I)=FMEAN0E(K,XCUT,XCUT_ERR,FSIGMA,YCUT_ERR(I))
                  ELSEIF(ICUT.EQ.2)THEN
                    YCUT(I)=
     +               FMEAN2E(K,XCUT,XCUT_ERR,3.0,FSIGMA,YCUT_ERR(I))
                  ELSEIF(ICUT.EQ.3)THEN
                    YCUT(I)=FMEDIAN1(K,XCUT)
                    YCUT_ERR(I)=0.0
                  END IF
                END DO
                IF(COPER.EQ.'+')THEN
                  IF(NBUFF0.LE.NMAXBUFF/2)THEN
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)+
     +                   YCUT_ERR(I)*YCUT_ERR(I))
                      END DO
                    END DO
                  END IF
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)+YCUT(I)
                    END DO
                  END DO
                ELSEIF(COPER.EQ.'-')THEN
                  IF(NBUFF0.LE.NMAXBUFF/2)THEN
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)+
     +                   YCUT_ERR(I)*YCUT_ERR(I))
                      END DO
                    END DO
                  END IF
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)-YCUT(I)
                    END DO
                  END DO
                ELSEIF(COPER.EQ.'*')THEN
                  IF(NBUFF0.LE.NMAXBUFF/2)THEN
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                   YCUT(I)*YCUT(I)+
     +                   IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                   YCUT_ERR(I)*YCUT_ERR(I))
                      END DO
                    END DO
                  END IF
                  DO I=NY1,NY2
                    DO J=NX1,NX2
                      IMAGEN(J,I,NBUFF0)=IMAGEN(J,I,NBUFF0)*YCUT(I)
                    END DO
                  END DO
                ELSEIF(COPER.EQ.'/')THEN
                  LDIVZERO=.FALSE.
                  DO I=NY1,NY2
                    IF(YCUT(I).EQ.0.0) LDIVZERO=.TRUE.
                  END DO
                  IF(LDIVZERO)THEN
                    WRITE(*,101) '***ERROR***'
                    WRITE(*,101) '=> Division by zero!'
                    CZERO(1:1)=READC('e[x]it or [c]ontinue','c','xc')
                    IF(CZERO.EQ.'c')THEN
                      FZERO=READF('Pixel value for infinity','0.0')
                      IF(NBUFF0.LE.NMAXBUFF/2)THEN
                        FZERO_ERR=READF(
     +                   'Pixel value for error at infinity','0.0')
                        DO I=NY1,NY2
                          IF(YCUT(I).EQ.0.0)THEN
                            DO J=NX1,NX2
                              IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=FZERO_ERR
                            END DO
                          ELSE
                            DO J=NX1,NX2
                              IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                         IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                         IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                         YCUT(I)*YCUT(I)+
     +                         IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                         YCUT_ERR(I)*YCUT_ERR(I))/
     +                         (YCUT(I)*YCUT(I))
                            END DO
                          END IF
                        END DO
                      END IF
                      DO I=NY1,NY2
                        IF(YCUT(I).EQ.0.0)THEN
                          DO J=NX1,NX2
                            IMAGEN(J,I,NBUFF0)=FZERO
                          END DO
                        ELSE
                          DO J=NX1,NX2
                            IMAGEN(J,I,NBUFF0)=
     +                       IMAGEN(J,I,NBUFF0)/YCUT(I)
                          END DO
                        END IF
                      END DO
                    END IF
                  ELSE
                    IF(NBUFF0.LE.NMAXBUFF/2)THEN
                      DO I=NY1,NY2
                        DO J=NX1,NX2
                          IMAGEN(J,I,NBUFF0+NMAXBUFF/2)=SQRT(
     +                     IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                     IMAGEN(J,I,NBUFF0+NMAXBUFF/2)*
     +                     YCUT(I)*YCUT(I)+
     +                     IMAGEN(J,I,NBUFF0)*IMAGEN(J,I,NBUFF0)*
     +                     YCUT_ERR(I)*YCUT_ERR(I))/
     +                     (YCUT(I)*YCUT(I))
                        END DO
                      END DO
                    END IF
                    DO I=NY1,NY2
                      DO J=NX1,NX2
                        IMAGEN(J,I,NBUFF0)=
     +                   IMAGEN(J,I,NBUFF0)/YCUT(I)
                      END DO
                    END DO
                  END IF
                END IF
              END IF
              LCANCEL=.TRUE.
            END IF
          END IF
C..............................................................................
        END DO
C------------------------------------------------------------------------------
C desactivamos todos los botones en la region de la calculadora
        DO J=1,4 !4 columnas de botones
          DO I=1,NMAXBUFF+2
            NB=10*(I-1)+NORIGEN+(J-1)
            CALL BUTTSEX(NB,.FALSE.)
          END DO
        END DO
        CALL RPGERASW(0.49,0.91,0.10,0.70,0)
        RETURN
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
