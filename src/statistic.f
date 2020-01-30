C
C******************************************************************************
C Calcula estadistica del rectangulo definido por NX1,NX2,NY1,NY2.
C Si LMEDIAN=.TRUE., tambien calcula la mediana (consume mucho tiempo)
C Si LSHOW=.TRUE., muestra los valores medidos
C Si LALL=.TRUE., usa toda la imagen
C Si LALL=.FALSE., usa solo puntos con sen~al entre BG y FG
C Si NBUFF=0, la estadistica se realiza sobre IMAGEN_(J,I) en lugar de 
C             sobre IMAGEN(J,I,NBUFF)
C Si LUSEMASK=.TRUE., entonces la estadistica se realiza sobre todos los pixels
C             seleccionados que no tengan como valor MASKVALUE
        SUBROUTINE STATISTIC(NBUFF,NX1,NX2,NY1,NY2,LMEDIAN,LSHOW,LALL,
     +   MASKVALUE,LUSEMASK)
        IMPLICIT NONE
        INTEGER NBUFF
        INTEGER NX1,NX2,NY1,NY2
        LOGICAL LMEDIAN,LSHOW,LALL
        REAL MASKVALUE
        LOGICAL LUSEMASK
C
        INCLUDE 'dimensions.inc'
C
        REAL FMEAN0
        REAL FMEAN0E
        REAL FMEAN2
        REAL FMEDIAN1
C
        INTEGER I,J,K
        INTEGER IFMIN,IFMAX,JFMIN,JFMAX
        INTEGER NPIX
        INTEGER IPIXEL(NXMAX*NYMAX), JPIXEL(NXMAX*NYMAX)
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL PIXEL(NXMAX*NYMAX),EPIXEL(NXMAX*NYMAX)
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        REAL EFMEAN
        REAL FMEANTSIGMA,FSIGMATSIGMA
        REAL Q25,Q75,FSIGMA_ROBUST
        REAL Q1587,Q8413
        REAL Q0250,Q9750
        REAL BG,FG
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/IMAGEN_
        COMMON/BLKESTADISTICA/FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        COMMON/BLKBGFG/BG,FG
C------------------------------------------------------------------------------
        EFMEAN=0.
        IF(LALL)THEN
          IF(NBUFF.EQ.0)THEN
            K=0
            DO I=NY1,NY2
              DO J=NX1,NX2
                IF(LUSEMASK)THEN
                  IF(IMAGEN_(J,I).NE.MASKVALUE)THEN
                    K=K+1
                    PIXEL(K)=IMAGEN_(J,I)
                    JPIXEL(K)=J
                    IPIXEL(K)=I
                  END IF
                ELSE
                  K=K+1
                  PIXEL(K)=IMAGEN_(J,I)
                  JPIXEL(K)=J
                  IPIXEL(K)=I
                END IF
              END DO
            END DO
          ELSE
            K=0
            DO I=NY1,NY2
              DO J=NX1,NX2
                IF(LUSEMASK)THEN
                  IF(IMAGEN(J,I,NBUFF).NE.MASKVALUE)THEN
                    K=K+1
                    PIXEL(K)=IMAGEN(J,I,NBUFF)
                    JPIXEL(K)=J
                    IPIXEL(K)=I
                  END IF
                ELSE
                  K=K+1
                  PIXEL(K)=IMAGEN(J,I,NBUFF)
                  JPIXEL(K)=J
                  IPIXEL(K)=I
                END IF
              END DO
            END DO
            IF(NBUFF.LE.NMAXBUFF/2)THEN
              K=0
              DO I=NY1,NY2
                DO J=NX1,NX2
                  IF(LUSEMASK)THEN
                    IF(IMAGEN(J,I,NBUFF).NE.MASKVALUE)THEN
                      K=K+1
                      EPIXEL(K)=IMAGEN(J,I,NBUFF+NMAXBUFF/2)
                    END IF
                  ELSE
                    K=K+1
                    EPIXEL(K)=IMAGEN(J,I,NBUFF+NMAXBUFF/2)
                  END IF
                END DO
              END DO
            END IF
          END IF
        ELSE
          IF(NBUFF.EQ.0)THEN
            K=0
            DO I=NY1,NY2
              DO J=NX1,NX2
                IF(LUSEMASK)THEN
                  IF(IMAGEN_(J,I).NE.MASKVALUE)THEN
                    IF((IMAGEN_(J,I).GE.BG).AND.
     +               (IMAGEN_(J,I).LE.FG))THEN
                      K=K+1
                      PIXEL(K)=IMAGEN_(J,I)
                      JPIXEL(K)=J
                      IPIXEL(K)=I
                    END IF
                  END IF
                ELSE
                  IF((IMAGEN_(J,I).GE.BG).AND.
     +             (IMAGEN_(J,I).LE.FG))THEN
                    K=K+1
                    PIXEL(K)=IMAGEN_(J,I)
                    JPIXEL(K)=J
                    IPIXEL(K)=I
                  END IF
                END IF
              END DO
            END DO
          ELSE
            K=0
            DO I=NY1,NY2
              DO J=NX1,NX2
                IF(LUSEMASK)THEN
                  IF(IMAGEN(J,I,NBUFF).NE.MASKVALUE)THEN
                    IF((IMAGEN(J,I,NBUFF).GE.BG).AND.
     +               (IMAGEN(J,I,NBUFF).LE.FG))THEN
                      K=K+1
                      PIXEL(K)=IMAGEN(J,I,NBUFF)
                      JPIXEL(K)=J
                      IPIXEL(K)=I
                    END IF
                  END IF
                ELSE
                  IF((IMAGEN(J,I,NBUFF).GE.BG).AND.
     +             (IMAGEN(J,I,NBUFF).LE.FG))THEN
                    K=K+1
                    PIXEL(K)=IMAGEN(J,I,NBUFF)
                    JPIXEL(K)=J
                    IPIXEL(K)=I
                  END IF
                END IF
              END DO
            END DO
            IF(NBUFF.LE.NMAXBUFF/2)THEN
              K=0
              DO I=NY1,NY2
                DO J=NX1,NX2
                  IF(LUSEMASK)THEN
                    IF(IMAGEN(J,I,NBUFF).NE.MASKVALUE)THEN
                      IF((IMAGEN(J,I,NBUFF).GE.BG).AND.
     +                 (IMAGEN(J,I,NBUFF).LE.FG))THEN
                        K=K+1
                        EPIXEL(K)=IMAGEN(J,I,NBUFF+NMAXBUFF/2)
                      END IF
                    END IF
                  ELSE
                    IF((IMAGEN(J,I,NBUFF).GE.BG).AND.
     +               (IMAGEN(J,I,NBUFF).LE.FG))THEN
                      K=K+1
                      EPIXEL(K)=IMAGEN(J,I,NBUFF+NMAXBUFF/2)
                    END IF
                  END IF
                END DO
              END DO
            END IF
          END IF
        END IF
        NPIX=K
C------------------------------------------------------------------------------
        IF(NPIX.EQ.0)THEN
          FMEAN=0.
          FSIGMA=0.
          FMEDIAN=0.
          FMIN=0.
          FMAX=0.
          WRITE(*,101) '***ERROR***'
          WRITE(*,100) '=> No. of pixels in STATISTICS=0.'
          WRITE(*,100) ' (press <CR>...)'
          READ(*,*)
        ELSEIF(NPIX.EQ.1)THEN
          IF(NBUFF.EQ.0)THEN
            FMEAN=IMAGEN_(NX1,NY1)
          ELSE
            FMEAN=IMAGEN(NX1,NY1,NBUFF)
          END IF
          FSIGMA=0.
          FMEDIAN=FMEAN
          FMIN=FMEAN
          FMAX=FMEAN
        ELSE
          IF(NBUFF.EQ.0)THEN
            FMEAN=FMEAN0(NPIX,PIXEL,FSIGMA)
          ELSE
            IF(NBUFF.LE.NMAXBUFF/2)THEN
              FMEAN=FMEAN0E(NPIX,PIXEL,EPIXEL,FSIGMA,EFMEAN)
            ELSE
              FMEAN=FMEAN0(NPIX,PIXEL,FSIGMA)
            END IF
          END IF
          FMIN=PIXEL(1)
          FMAX=FMIN
          JFMIN=1
          IFMIN=1
          JFMAX=1
          IFMAX=1
          DO K=2,NPIX
            IF(PIXEL(K).LT.FMIN)THEN
              FMIN=PIXEL(K)
              JFMIN=JPIXEL(K)
              IFMIN=IPIXEL(K)
            END IF
            IF(PIXEL(K).GT.FMAX)THEN
              FMAX=PIXEL(K)
              JFMAX=JPIXEL(K)
              IFMAX=IPIXEL(K)
            END IF
          END DO
          IF(LSHOW)THEN
            WRITE(*,*)
            WRITE(*,100) '=> Number of pixels..........: '
            WRITE(*,*) NPIX
            WRITE(*,100) '=> Minimum, at (X,Y).........: '
            WRITE(*,*) FMIN,JFMIN,IFMIN
            WRITE(*,100) '=> Maximum, at (X,Y).........: '
            WRITE(*,*) FMAX,JFMAX,IFMAX
            WRITE(*,100) '=> Mean & error..............: '
            WRITE(*,*) FMEAN,EFMEAN
            WRITE(*,100) '=> Total no. of counts & err.: '
            WRITE(*,*) FMEAN*REAL(NPIX),EFMEAN*REAL(NPIX)
            WRITE(*,100) '=> StDev.....................: '
            WRITE(*,*) FSIGMA
          END IF
          IF(LMEDIAN)THEN
            FMEDIAN=FMEDIAN1(NPIX,PIXEL)
            IF(LSHOW)THEN
              Q25=PIXEL(INT(REAL(NPIX)*0.25+0.5))
              Q75=PIXEL(INT(REAL(NPIX)*0.75+0.5))
              FSIGMA_ROBUST=0.7413*(Q75-Q25)
              WRITE(*,100) '=> StDev (robust)............: '
              WRITE(*,*) FSIGMA_ROBUST
              WRITE(*,100) '=> Median....................: '
              WRITE(*,*) FMEDIAN
              Q0250=PIXEL(INT(REAL(NPIX)*0.0250+0.5))
              Q1587=PIXEL(INT(REAL(NPIX)*0.1587+0.5))
              Q8413=PIXEL(INT(REAL(NPIX)*0.8413+0.5))
              Q9750=PIXEL(INT(REAL(NPIX)*0.9750+0.5))
              WRITE(*,100) '=> Q 02.50%..................: '
              WRITE(*,*) Q0250
              WRITE(*,100) '=> Q 15.87%..................: '
              WRITE(*,*) Q1587
              WRITE(*,100) '=> Q 84.13%..................: '
              WRITE(*,*) Q8413
              WRITE(*,100) '=> Q 97.50%..................: '
              WRITE(*,*) Q9750
            END IF
            FMEANTSIGMA=FMEAN2(NPIX,PIXEL,3.0,FSIGMATSIGMA)   
            IF(LSHOW)THEN
              WRITE(*,100) '=> Mean  (removing > 3 sigma): '
              WRITE(*,*) FMEANTSIGMA
              WRITE(*,100) '=> StDev (removing > 3 sigma): '
              WRITE(*,*) FSIGMATSIGMA
            END IF
          ELSE
            FMEDIAN=0.
          END IF
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
