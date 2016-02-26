C
C******************************************************************************
C Calcula estadistica del rectangulo definido por NX1,NX2,NY1,NY2 en un box-9.
C Si NBUFF=0, la estadistica se realiza sobre IMAGEN_(J,I) en lugar de 
C             sobre IMAGEN(J,I,NBUFF)
        SUBROUTINE STATISTICB9(NFRAMES_,NBUFF,FMEAN,FSIGMA)
        IMPLICIT NONE
        INTEGER NFRAMES_
        INTEGER NBUFF
        REAL FMEAN,FSIGMA
C
        INCLUDE 'dimensions.inc'
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        REAL FMEAN0
C
        INTEGER I,J,K
        INTEGER NPIX,NF
        INTEGER DI(NBOXMAX),DJ(NBOXMAX)
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF),PIXEL(NXMAX*NYMAX)
        REAL IMAGEN_(NXMAX,NYMAX)
C
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKIMAGEN1_/IMAGEN_
C------------------------------------------------------------------------------
C Note: the pattern of the frames in box-9 is the following:
C       6 9 4
C       3 1 7
C       8 5 2
C offsets of each frame in the 768x768 composite mosaic
        DATA (DI(NF),DJ(NF),NF=1,NBOXMAX) /
     +   256,256,  !1
     +   000,512,  !2
     +   256,000,  !3
     +   512,512,  !4
     +   000,256,  !5
     +   512,000,  !6
     +   256,512,  !7
     +   000,000,  !8
     +   512,256/  !9
C------------------------------------------------------------------------------
        K=0
        DO NF=1,NFRAMES_
          IF(NBUFF.EQ.0)THEN
            DO I=1,256
              DO J=1,256
                K=K+1
                PIXEL(K)=IMAGEN_(J+DJ(NF),I+DI(NF))
              END DO
            END DO
          ELSE
            DO I=1,256
              DO J=1,256
                K=K+1
                PIXEL(K)=IMAGEN(J+DJ(NF),I+DI(NF),NBUFF)
              END DO
            END DO
          END IF
        END DO
        NPIX=K
C------------------------------------------------------------------------------
        FMEAN=FMEAN0(NPIX,PIXEL,FSIGMA)
C------------------------------------------------------------------------------
        END
