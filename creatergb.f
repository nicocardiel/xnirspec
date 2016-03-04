        SUBROUTINE CREATERGB(NCBUFF)
        IMPLICIT NONE
        INTEGER NCBUFF
C
        INCLUDE 'dimensions.inc'
C
        INTEGER SYSTEMFUNCTION
        INTEGER READILIM
        INTEGER TRUEBEG,TRUELEN
        REAL READF
        CHARACTER*255 READC
C
        INTEGER I,J,K,L,L1,L2
        INTEGER NX1,NX2,NY1,NY2
        INTEGER NBUFF_R,NBUFF_G,NBUFF_B
        INTEGER NX,NY
        INTEGER IPIXELR,IPIXELG,IPIXELB
        INTEGER ISYSTEM
        INTEGER NEXTINFO
        INTEGER NCOLOROUT
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL BG,FG
        REAL PIXELR,PIXELG,PIXELB
        REAL FRGB_R,FRGB_G,FRGB_B
        REAL FLOGRGB
        REAL PIXELMAX
        CHARACTER*1 CRGBSATUR,CENHANCE_RB
        CHARACTER*1 COVER
        CHARACTER*1 CLOGRGB
        CHARACTER*50 CDUMMY
        CHARACTER*80 PPMFILE
        LOGICAL LOGFILE
        LOGICAL LOUTR,LOUTG,LOUTB
C
        COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKBGFG/BG,FG
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKRGBBUFF/NBUFF_R,NBUFF_G,NBUFF_B
        COMMON/BLKRGBFACT/FRGB_R,FRGB_G,FRGB_B,FLOGRGB
        COMMON/BLKRGBSATUR/CRGBSATUR,CENHANCE_RB,CLOGRGB
C------------------------------------------------------------------------------
        LOGFILE=.TRUE.
        L1=0 !evita un warning the compilacion
        L2=0 !evita un warning the compilacion
        DO WHILE(LOGFILE)
          PPMFILE(1:80)=READC('Output ppm file name','pgplot.ppm','@')
          L1=TRUEBEG(PPMFILE)
          L2=TRUELEN(PPMFILE)
          INQUIRE(FILE=PPMFILE(L1:L2),EXIST=LOGFILE)
          IF(LOGFILE)THEN
            IF(PPMFILE(L1:L2).EQ.'pgplot.ppm')THEN
              COVER='y' !no pedimos confirmacion en este caso
            ELSE
              WRITE(*,101) 'WARNING: this file already exist.'
              COVER(1:1)=
     +         READC('Do you want to overwrite it (y/n)','y','yn')
            END IF
            IF(COVER.EQ.'y')THEN
              ISYSTEM=SYSTEMFUNCTION('rm '//PPMFILE(L1:L2))
              LOGFILE=.FALSE.
            END IF
          END IF
        END DO
C
        IF(NBUFF_R.EQ.-1) NBUFF_R=NCBUFF
        WRITE(CDUMMY,*) NBUFF_R
        NBUFF_R=READILIM('Buffer # for RED',CDUMMY,1,NMAXBUFF)
        IF(NBUFF_G.EQ.-1) NBUFF_G=NBUFF_R+1
        WRITE(CDUMMY,*) NBUFF_G
        NBUFF_G=READILIM('Buffer # for GREEN (0=(red+blue)/2)',
     +   CDUMMY,0,NMAXBUFF)
        IF(NBUFF_B.EQ.-1)THEN
          IF(NBUFF_G.EQ.0)THEN
            NBUFF_B=NBUFF_R+1
          ELSE
            NBUFF_B=NBUFF_G+1
          END IF
        END IF
        WRITE(CDUMMY,*) NBUFF_B
        NBUFF_B=READILIM('Buffer # for BLUE',CDUMMY,1,NMAXBUFF)
C
        WRITE(CDUMMY,*) FRGB_R
        FRGB_R=READF('Factor to scale RED image',CDUMMY)
        IF(NBUFF_G.NE.0)THEN
          WRITE(CDUMMY,*) FRGB_G
          FRGB_G=READF('Factor to scale GREEN image',CDUMMY)
        END IF
        WRITE(CDUMMY,*) FRGB_B
        FRGB_B=READF('Factor to scale BLUE image',CDUMMY)
C
        CRGBSATUR(1:1)=
     +   READC('Are you using saturation correction',
     +   CRGBSATUR,'yn')
C
        NCOLOROUT=READILIM('Color number for region outside image '//
     +   '(gray scale)','0',0,255)
C
        CENHANCE_RB(1:1)=READC('Enhance Blue and Red',CENHANCE_RB,'yn')
C
        CLOGRGB(1:1)=READC('Use a logarithmic scale',CLOGRGB,'yn')
        IF(CLOGRGB.EQ.'y')THEN
          WRITE(CDUMMY,*) FLOGRGB
          FLOGRGB=READF('Factor for logarithmic scale',CDUMMY)
        END IF
C------------------------------------------------------------------------------
C Abrimos fichero y escribimos cabecera
        OPEN(20,FILE=PPMFILE,STATUS='UNKNOWN',FORM='FORMATTED')
        WRITE(20,101) 'P3'
        WRITE(20,101) '# CREATOR: xnirspec Version 2.0'
        NX=NX2-NX1+1
        NY=NY2-NY1+1
        WRITE(CDUMMY,*) NX
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(20,100) CDUMMY(1:L)
        WRITE(CDUMMY,*) NY
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(20,101) ' '//CDUMMY(1:L)
        WRITE(20,101) '255'
C escribimos imagen
        K=0
        NEXTINFO=0
        WRITE(*,101) 'Creating ppm file...'
        DO I=NY2,NY1,-1 !OJO: hay que escribir en orden inverso en Y
          DO J=NX1,NX2
            K=K+1
            LOUTR=.FALSE.
            LOUTG=.FALSE.
            LOUTB=.FALSE.
            IF(NBUFF_R.GT.NMAXBUFF/2)THEN
              PIXELR=IMAGEN(J,I,NBUFF_R)*FRGB_R
            ELSE
              IF(IMAGEN(J,I,NBUFF_R+NMAXBUFF/2).GE.-0.5)THEN
                PIXELR=IMAGEN(J,I,NBUFF_R)*FRGB_R
              ELSE
                PIXELR=BG
                LOUTR=.TRUE.
              END IF
            END IF
            IF(NBUFF_B.GT.NMAXBUFF/2)THEN
              PIXELB=IMAGEN(J,I,NBUFF_B)*FRGB_B
            ELSE
              IF(IMAGEN(J,I,NBUFF_B+NMAXBUFF/2).GE.-0.5)THEN
                PIXELB=IMAGEN(J,I,NBUFF_B)*FRGB_B
              ELSE
                PIXELB=BG
                LOUTB=.TRUE.
              END IF
            END IF
            IF(NBUFF_G.EQ.0)THEN
              PIXELG=(PIXELR+PIXELB)/2.
            ELSE
              IF(NBUFF_G.GT.NMAXBUFF/2)THEN
                PIXELG=IMAGEN(J,I,NBUFF_G)*FRGB_G
              ELSE
                IF(IMAGEN(J,I,NBUFF_G+NMAXBUFF/2).GE.-0.5)THEN
                  PIXELG=IMAGEN(J,I,NBUFF_G)*FRGB_G
                ELSE
                  PIXELG=BG
                  LOUTG=.TRUE.
                END IF
              END IF
            END IF
            IF(NBUFF_G.EQ.0)THEN
              IF(LOUTR.AND.LOUTB) LOUTG=.TRUE.
            END IF
            IF(LOUTR.AND.LOUTG.AND.LOUTB)THEN
              IPIXELR=NCOLOROUT
              IPIXELG=NCOLOROUT
              IPIXELB=NCOLOROUT
            ELSE
              IF(CENHANCE_RB.EQ.'y')THEN
                IF(PIXELR.GT.PIXELB)THEN
                  PIXELR=PIXELR+PIXELB
                  PIXELB=0.
                ELSEIF(PIXELR.LT.PIXELB)THEN
                  PIXELB=PIXELB+PIXELR
                  PIXELR=0.
                END IF
              END IF
              IF(CLOGRGB.EQ.'y')THEN
                IF(PIXELR.GT.BG)THEN
                  PIXELR=(ALOG10((PIXELR-BG)/(FG-BG))+FLOGRGB)/
     +             FLOGRGB*255.
                ELSE
                  PIXELR=0.
                END IF
                IF(PIXELG.GT.BG)THEN
                  PIXELG=(ALOG10((PIXELG-BG)/(FG-BG))+FLOGRGB)/
     +             FLOGRGB*255.
                ELSE
                  PIXELG=0.
                END IF
                IF(PIXELB.GT.BG)THEN
                  PIXELB=(ALOG10((PIXELB-BG)/(FG-BG))+FLOGRGB)/
     +             FLOGRGB*255.
                ELSE
                  PIXELB=0.
                END IF
              ELSE
                PIXELR=(PIXELR-BG)/(FG-BG)*255.
                PIXELG=(PIXELG-BG)/(FG-BG)*255.
                PIXELB=(PIXELB-BG)/(FG-BG)*255.
              END IF
              IF(CRGBSATUR.EQ.'y')THEN
                PIXELMAX=PIXELR
                IF(PIXELMAX.LT.PIXELG) PIXELMAX=PIXELG
                IF(PIXELMAX.LT.PIXELB) PIXELMAX=PIXELB
                IF(PIXELMAX.GT.255.)THEN !re-escalamos en caso de saturacion
                  PIXELR=PIXELR*255./PIXELMAX
                  PIXELG=PIXELG*255./PIXELMAX
                  PIXELB=PIXELB*255./PIXELMAX
                END IF
              END IF
              IPIXELR=NINT(PIXELR)
              IF(IPIXELR.LT.0) IPIXELR=0
              IF(IPIXELR.GT.255) IPIXELR=255
              IPIXELG=NINT(PIXELG)
              IF(IPIXELG.LT.0) IPIXELG=0
              IF(IPIXELG.GT.255) IPIXELG=255
              IPIXELB=NINT(PIXELB)
              IF(IPIXELB.LT.0) IPIXELB=0
              IF(IPIXELB.GT.255) IPIXELB=255
            END IF
            IF(K.EQ.1)THEN
              WRITE(20,'(I3,1X,I3,1X,I3,$)') IPIXELR,IPIXELG,IPIXELB
            ELSE
              WRITE(20,'(1X,I3,1X,I3,1X,I3,$)') IPIXELR,IPIXELG,IPIXELB
            END IF
            IF(K.EQ.5)THEN
              K=0
              WRITE(20,101)
            END IF
          END DO
          CALL SHOWPERC(1,NY2-NY1+1,1,NY2-I+1,NEXTINFO)
        END DO
C en caso necesario, completamos ultima linea
        IF(MOD(NX*NY,5).NE.0)THEN
          DO I=1,5-MOD(NX*NY,5)
            WRITE(20,'(1X,I3,1X,I3,1X,I3,$)') 0,0,0
          END DO
          WRITE(20,101)
        END IF
C cerramos fichero
        CLOSE(20)
C abrimos imagen con xv
!       WRITE(*,100) 'Launching viewer...'
!       ISYSTEM=SYSTEMFUNCTION('xv '//PPMFILE(L1:L2)//' &')
!       WRITE(*,*)
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
