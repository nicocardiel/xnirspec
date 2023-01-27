! Desplaza la imagen NBUFF2 para que coincida con la imagen NBUFF1 usando
! el offset que se mide al usar un mismo objeto en las dos imágenes
!
        SUBROUTINE SHIFTIMAGE
        IMPLICIT NONE
!
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER READILIM
!
        INTEGER NCBUFF
        INTEGER NBUFF1,NBUFF2
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER NAXIS1,NAXIS2
        INTEGER NAXIS1_,NAXIS2_
        INTEGER I,J,I1,J1,I2,J2
        INTEGER I0,J0
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)                        !current image
        REAL IMAGEN_(NXMAX,NYMAX)
        REAL XC,YC                                              !mouse location
        REAL X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE     !coefficients of centroid fit
        REAL EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE !errors
        REAL X0_,Y0_
        REAL XOFFSET,YOFFSET
        REAL FJ0,FI0
        CHARACTER*1 CH
        CHARACTER*50 CDUMMY
        LOGICAL LOOP
! common blocks
        COMMON/BLKIMAGEN2/NCBUFF                      !numero del buffer actual
        COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKIMAGEN1_/IMAGEN_              !es global para ahorrar memoria
        COMMON/BLKNAXIS/NAXIS
!------------------------------------------------------------------------------
        I2=0 !avoid compilation warning
        J2=0 !avoid compilation warning
!------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,101) '============================================'
        WRITE(*,101) 'Note: this option overwrites the image to be'
        WRITE(*,101) '      shifted. Error frames are NOT handled;'
        WRITE(*,101) '      if you need this use the shift option'
        WRITE(*,101) '      in calculator.'
        WRITE(*,101) '============================================'
! pedimos buffer de referencia
        NBUFF1=NCBUFF      !por defecto suponemos que va a ser el buffer activo
        LOOP=.TRUE.
        DO WHILE(LOOP)
          WRITE(CDUMMY,*) NBUFF1
          NBUFF1=READILIM('Buffer # of REFERENCE image (0=exit)',CDUMMY,0,NMAXBUFF/2)
          IF(NBUFF1.EQ.0) RETURN
          NAXIS1 = NAXIS(1,NBUFF1)
          NAXIS2 = NAXIS(2,NBUFF1)
          IF((NAXIS1.EQ.0).AND.(NAXIS2.EQ.0))THEN
            WRITE(*,101) 'Invalid buffer (undefined!). Try again.'
          ELSE
            LOOP=.FALSE.
          END IF
        END DO
! pedimos buffer con imagen a ser desplazada
        LOOP=.TRUE.
        DO WHILE(LOOP)
          NBUFF2=READILIM('Buffer # of image to be SHIFTED (0=exit)','@',0,NMAXBUFF/2)
          IF(NBUFF2.EQ.0) RETURN
          IF(NBUFF2.EQ.NBUFF1)THEN
            WRITE(*,100) 'ERROR: you must select a different buffer!'
            WRITE(*,101) ' Try again.'
          ELSE
            NAXIS1_ = NAXIS(1,NBUFF2)
            NAXIS2_ = NAXIS(2,NBUFF2)
            IF((NAXIS1.NE.NAXIS1_).OR.(NAXIS2.NE.NAXIS2_))THEN
              WRITE(*,100) '>>> Buffer#, naxis1, naxis2: '
              WRITE(*,*) NBUFF1, NAXIS1, NAXIS2
              WRITE(*,100) '>>> Buffer#, naxis1, naxis2: '
              WRITE(*,*) NBUFF2, NAXIS1_, NAXIS2_
              WRITE(*,101) 'Buffer dimensions do not match! Try again.'
            ELSE
              LOOP=.FALSE.
            END IF
          END IF
        END DO
!------------------------------------------------------------------------------ 
! dibujamos primera imagen y medimos objeto
        CALL SUBLOOK(.TRUE.,NBUFF1,.FALSE.)
        WRITE(*,*)
        WRITE(*,101) 'Select object in REFERENCE image'
        CALL CENTROID(NBUFF1,CH,XC,YC,.FALSE.,X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE,EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,.TRUE.)
! dibujamos segunda imagen y medimos el mismo objeto
        CALL SUBLOOK(.TRUE.,NBUFF2,.FALSE.)
        WRITE(*,*)
        WRITE(*,101) 'Select object in image to be SHIFTED'
        CALL CENTROID(NBUFF2,CH,XC,YC,.FALSE.,X0_,Y0_,SIGMAX,SIGMAY,BETA,AMP,CTE,EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,.TRUE.)
! calculamos offset
        XOFFSET=X0-X0_
        YOFFSET=Y0-Y0_
        WRITE(*,*)
        WRITE(*,100) '>>> X offset: '
        WRITE(*,*) XOFFSET
        WRITE(*,100) '>>> Y offset: '
        WRITE(*,*) YOFFSET
        WRITE(*,100) 'Applying previous offset'
!------------------------------------------------------------------------------ 
! desplazamos la imagen (el codigo aqui lo copio de la subrutina SUBIMATH,
! aunque no manejamos imágenes de errores y la imagen de destino se
! sobreescribe; es decir, es una version simplificada)
        J0=INT(XOFFSET)
        FJ0=XOFFSET-REAL(J0)
        I0=INT(YOFFSET)
        FI0=YOFFSET-REAL(I0)
! inicializamos imagen de destino a cero
        DO I=1,NAXIS2
          DO J=1,NAXIS1
            IMAGEN_(J,I)=0.
          END DO
        END DO
! aplicamos el offset a la imagen de datos (I,J son el pixel en la imagen
! inicial, mientras que I1,J y I2,J son los dos posibles pixels a los que
! movemos la sen~al en la nueva imagen)
        DO I=1,NAXIS2
          I1=I+I0
          IF(FI0.EQ.0.0)THEN
            I2=0
          ELSEIF(FI0.GT.0.0)THEN
            I2=I1+1
          ELSEIF(FI0.LT.0.0)THEN
            I2=I1-1
          END IF
          DO J=1,NAXIS1
            J1=J+J0
            IF(FJ0.EQ.0.0)THEN
              J2=0
            ELSEIF(FJ0.GT.0.0)THEN
              J2=J1+1
            ELSEIF(FJ0.LT.0.0)THEN
              J2=J1-1
            END IF
            IF((J1.GE.1).AND.(J1.LE.NAXIS1))THEN
              IF((I1.GE.1).AND.(I1.LE.NAXIS2))THEN
                IMAGEN_(J1,I1)=IMAGEN_(J1,I1)+IMAGEN(J,I,NBUFF2)*(1.-ABS(FI0))*(1.-ABS(FJ0))
              END IF
              IF((I2.GE.1).AND.(I2.LE.NAXIS2))THEN
                IMAGEN_(J1,I2)=IMAGEN_(J1,I2)+IMAGEN(J,I,NBUFF2)*ABS(FI0)*(1.-ABS(FJ0))
              END IF
            END IF
            IF((J2.GE.1).AND.(J2.LE.NAXIS1))THEN
              IF((I1.GE.1).AND.(I1.LE.NAXIS2))THEN
                IMAGEN_(J2,I1)=IMAGEN_(J2,I1)+IMAGEN(J,I,NBUFF2)*(1.-ABS(FI0))*ABS(FJ0)
              END IF
              IF((I2.GE.1).AND.(I2.LE.NAXIS2))THEN
                IMAGEN_(J2,I2)=IMAGEN_(J2,I2)+IMAGEN(J,I,NBUFF2)*ABS(FI0)*ABS(FJ0)
              END IF
            END IF
          END DO
        END DO
! introducimos en el buffer de destino el resultado
        DO I=1,NAXIS2
          DO J=1,NAXIS1
            IMAGEN(J,I,NBUFF2)=IMAGEN_(J,I)
          END DO
        END DO
! dibujamos la segunda imagen en su nueva posición
        CALL SUBLOOK(.TRUE.,NBUFF2,.FALSE.)
!------------------------------------------------------------------------------ 
100     FORMAT(A,$)
101     FORMAT(A)
        END
