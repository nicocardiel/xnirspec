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
        CHARACTER*255 INFILE(NMAXBUFF)
        LOGICAL LDS9REG(NMAXBUFF)
        LOGICAL LASCREG(NMAXBUFF)
! common blocks
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
        COMMON/BLKBGFG/BG,FG
        COMMON/BLKJUST/JUST
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKINFILE/INFILE
        COMMON/BLKLDS9REG/LDS9REG
        COMMON/BLKLASCREG/LASCREG
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
            CALL RPGERASW(0.4,1.0,0.0,0.8,0)
            CALL RPGENV(REAL(NX1)-0.6,REAL(NX2)+0.6,REAL(NY1)-0.6,REAL(NY2)+0.6,JUST,-2)
          END IF
        END IF
        IF(LREVERSE)THEN
          CALL PALETTE(1)
          CALL PGGRAY(IMAGEN(1,1,NCBUFF),NXMAX,NYMAX,NX1,NX2,NY1,NY2,BG,FG,TR)
          CALL PALETTE(3)
        ELSE
          CALL PGIMAG(IMAGEN(1,1,NCBUFF),NXMAX,NYMAX,NX1,NX2,NY1,NY2,FG,BG,TR)
        END IF
        CALL PGBOX('IBCTSN',0.0,0,'IBCTSN',0.0,0)
        IF(.NOT.LOVER)THEN
          CALL PGSCI(5)
          CALL PGMTXT('T',1.0,0.5,0.5,INFILE(NCBUFF))
          CALL PGSCI(1)
        END IF
!------------------------------------------------------------------------------
        IF(LDS9REG(NCBUFF)) CALL PLOTDS9REG(NCBUFF)
        IF(LASCREG(NCBUFF)) CALL PLOTASCREG(NCBUFF)
!------------------------------------------------------------------------------
        END
