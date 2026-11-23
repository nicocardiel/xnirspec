! if LOVER=.TRUE. then overplot
! if LREVERSE=.TRUE. reverse BG and FG (useful for postscript plots)
        SUBROUTINE SUBLOOK(LOVER,NCBUFF,LREVERSE)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
! subroutine arguments
        LOGICAL LOVER
        INTEGER NCBUFF
        LOGICAL LREVERSE
!
        INCLUDE 'dimensions.inc'
!
        INTEGER JUST
        INTEGER NX1,NX2,NY1,NY2             !limits of current displayed region
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL BG,FG
        REAL TR(6)                          !auxiliary matrix for PGIMAG/PGGRAY
        REAL OLD_CH                                !previous character height
        REAL XV1,XV2,YV1,YV2                          !viewport (NDC)
        INTEGER NPALETTE              !current color palette (see PALETTE)
        CHARACTER*255 INFILE(NMAXBUFF)
        LOGICAL LDS9REG(NMAXBUFF)
        LOGICAL LASCREG(NMAXBUFF)
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)
        LOGICAL LWAVECAL(NMAXBUFF)
! common blocks
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKBGFG/BG,FG
        COMMON/BLKJUST/JUST
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKINFILE/INFILE
        COMMON/BLKLDS9REG/LDS9REG
        COMMON/BLKLASCREG/LASCREG
        COMMON/BLKPALETTE/NPALETTE
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1         !wavelength calibration
        COMMON/BLKWAVECAL2/LWAVECAL
!------------------------------------------------------------------------------
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
!------------------------------------------------------------------------------
        IF(.NOT.LOVER)THEN
          IF(LREVERSE)THEN !postscript file
            CALL PGENV(REAL(NX1)-0.6,REAL(NX2)+0.6,REAL(NY1)-0.6,REAL(NY2)+0.6,JUST,-2)
          ELSE
            CALL RPGERASW(0.4,1.0,0.0,0.76,0)
            CALL RPGENV(REAL(NX1)-0.6,REAL(NX2)+0.6,REAL(NY1)-0.6,REAL(NY2)+0.6,JUST,-2)
          END IF
        ELSEIF(.NOT.LREVERSE)THEN
! the image is redrawn over the previous one (e.g. another buffer of the
! same size): erase the strip between the top edge of the image and the
! row of buttons, which holds the title or the wavelength scale of the
! previous image (both are drawn again below)
          CALL PGQVP(0,XV1,XV2,YV1,YV2)
          CALL RPGERASW(0.4,1.0,YV2,0.76,0)
        END IF
        IF(LREVERSE)THEN
          CALL PALETTE(4)                                   !gray scale
          CALL PGGRAY(IMAGEN(1,1,NCBUFF),NXMAX,NYMAX,NX1,NX2,NY1,NY2,BG,FG,TR)
          CALL PALETTE(NPALETTE)        !restore the palette chosen by the user
        ELSE
          CALL PGIMAG(IMAGEN(1,1,NCBUFF),NXMAX,NYMAX,NX1,NX2,NY1,NY2,FG,BG,TR)
        END IF
        IF(LWAVECAL(NCBUFF))THEN
! wavelength calibrated image: wavelength scale along the top edge
          CALL PGBOX('IBTSN',0.0,0,'IBCTSN',0.0,0)
          CALL WVSCALEBOX(CRPIX1(NCBUFF),CRVAL1(NCBUFF),CDELT1(NCBUFF))
        ELSE
          CALL PGBOX('IBCTSN',0.0,0,'IBCTSN',0.0,0)
        END IF
        IF((.NOT.LOVER).OR.(.NOT.LREVERSE))THEN
! title with the same character height as in the plot of cuts (SUBPLOT), so
! that it does not reach the row of buttons above the image; with the
! wavelength scale there is no room above the image, and the title is
! written inside it (upper-left corner). The title is also written when the
! image is redrawn (LOVER=.TRUE.), since the strip above the image has been
! erased (or the image has been drawn over the title)
          CALL PGQCH(OLD_CH)
          CALL PGSCH(0.8)
          CALL PGSCI(5)
          IF(LWAVECAL(NCBUFF))THEN
            CALL PGMTXT('T',-1.5,0.02,0.0,INFILE(NCBUFF))
          ELSE
            CALL PGMTXT('T',1.0,0.5,0.5,INFILE(NCBUFF))
          END IF
          CALL PGSCI(1)
          CALL PGSCH(OLD_CH)
        END IF
!------------------------------------------------------------------------------
        IF(LDS9REG(NCBUFF)) CALL PLOTDS9REG(NCBUFF)
        IF(LASCREG(NCBUFF)) CALL PLOTASCREG(NCBUFF)
!------------------------------------------------------------------------------
        END
