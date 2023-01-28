!******************************************************************************
! Initial version 26-March-2002        (c) Nicolas Cardiel: cardiel@ucolick.org
!                                                         ncl@astrax.fis.ucm.es
!******************************************************************************
! g77 -o xnirspec *.o -L$PGPLOT_DIR -L$FITSIO_DIR -L/usr/X11/lib \
! -lpgplot -lcfitsio -lX11
!
        PROGRAM XNIRSPEC
        USE Dynamic_Array_IMAGEN
        USE Dynamic_Array_IMAGEN_
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
        INCLUDE 'interface_imagen_.inc'
! parameters
        INCLUDE 'dimensions.inc'
        INCLUDE 'largest.inc'
!
        INTEGER NMAPMAX
        PARAMETER (NMAPMAX=3) !second-degree bivariate approximation
        INTEGER NECUAMAX
        PARAMETER (NECUAMAX=(NMAPMAX+1)*NMAPMAX)
        INTEGER NWIDTHMAX
        PARAMETER (NWIDTHMAX=101)
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
! functions
        INTEGER TRUEBEG
        INTEGER TRUELEN
        INTEGER SYSTEMFUNCTION
        INTEGER READILIM
        INTEGER READI
        REAL READF
        REAL FPOLY
        CHARACTER*255 READC
! variables que definen el boundary
! NOTA: Por razones historicas, los polinomios que definen el boundary se 
! denominan BL y BA por Boundary Lines, y Boundary Arc lines. 
! Cuando se decidio rotar las imagenes para tener dispersion en el 
! eje X y direccion espacial en el eje Y, BL se convirtio en polinomios 
! en la direccion espectral de la forma y=f(x), y BA en polinomios en la 
! direccion espacial de la forma x=g(y).
        INTEGER NXYMAX_
        INTEGER NDEGBL(NXYMAX),NDEGBA(NXYMAX)  !polynomial degrees for boundary
        INTEGER NDEGBL00,NDEGBA00                 !maximum of NDEGBL and NDEGBA
        INTEGER NLINBL,NLINBA            !number of BL and BA lines in boundary
        REAL COEFFBL(20,NXYMAX),COEFFBA(20,NXYMAX)      !pol. coef. in boundary
        REAL XMINBL(NXYMAX),XMAXBL(NXYMAX)
        REAL YMINBL(NXYMAX),YMAXBL(NXYMAX)
        REAL XMINBA(NXYMAX),XMAXBA(NXYMAX)
        REAL YMINBA(NXYMAX),YMAXBA(NXYMAX)
        LOGICAL LBOUNDARY                  !if .TRUE. boundary has been defined
        LOGICAL LBOUNDL,LBOUNDA          !.TRUE. if BL and BA lines are defined
        LOGICAL LPBOUND                                !if .TRUE. plot boundary
        LOGICAL LMAPPING                       !if .TRUE., mapping is available
        LOGICAL LPMAPP                             !if .TRUE. plot mapping grid
        LOGICAL LPWAVE                               !if .TRUE. plot wavelength
! otras variables
        INTEGER I1,I2,J1,J2                         !borders of rectified image
        INTEGER NCBUFF,NCBUFF_                           !current buffer in use
        INTEGER NEWBUFF,NEWBUFF_                         !new buffer to be used
        INTEGER NBUFFBOX9
        INTEGER NY_CUT         !corte inicial de la primera imagen que se carga
        INTEGER I,J,K,L,II,JJ                                 !generic counters
        INTEGER IBOX9                         !generic counter for loop in BOX9
        INTEGER PGOPEN                                         !PGPLOT function
        INTEGER ID,IDNEW              !identification number of graphic display
        INTEGER NB,NBLOCAL                           !number of selected button
        INTEGER L1,L2,L3,L4              !TRUEBEG and TRUELEN of generic string
        INTEGER ISYSTEM              !returned value of function SYSTEMFUNCTION
        INTEGER NAXIS(2,NMAXBUFF)                             !image dimensions
        INTEGER NAXIS1_,NAXIS2_                               !image dimensions
        INTEGER DAXIS1_,DAXIS2_                        !offsets to resize image
        INTEGER NAXISFRAME(2,9,NMAXBUFF)     !dimensions of each frame in box-9
        INTEGER DI(9),DJ(9)           !offsets in box-9 mosaic (in frame units)
        INTEGER NFRAMES(NMAXBUFF)                  !real number of frames/image
        INTEGER NSIZEFB9(NMAXBUFF)       !maximum size of frame in box-9 mosaic
        INTEGER JUST                                      !JUST value in PGPLOT
        INTEGER NX1,NX2,NY1,NY2             !limits of current displayed region
        INTEGER NX1_,NX2_,NY1_,NY2_
        INTEGER IXC,IYC                     !current mouse location (in pixels)
        INTEGER NCOLOR1,NCOLOR2,NCOLOR3,NCOLOR4,NCOLOR5,NCOLOR6  !button colors
        INTEGER MODECUT                          !1, 2, 3 = single, region, all
        INTEGER K1,K2                              !summed region in X or Y cut
        INTEGER NSTACK                           !no. of points stored in stack
        INTEGER NDEG                                         !polynomial degree
        INTEGER NDEGWV             !polynomial degree of wavelength calibration
        INTEGER NMAP                            !bivariate approximation degree
        INTEGER IROTATE                                         !rotation angle
        INTEGER SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID   !extremos del mapping
        INTEGER SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG !extrapola. del mapping
        INTEGER NBOX                                   !box size for statistics
        INTEGER NNX1,NNX2,NNY1,NNY2                 !box corners for statistics
        INTEGER NNX1_,NNY1_         !initial box corner for statistics in box 9
        INTEGER NYFINDAL                     !Y grid looking for spectral lines
        INTEGER NINLINE                !no. of peaks/line found (max. NYFINDAL)
        INTEGER NWIDTH            !no. pixels to search for peaks (odd, min. 3)
        INTEGER NGRID                      !grid width (pixels) to plot mapping
        INTEGER IRESAMPL          !resampling method (1=approximate,2=accurate)
        INTEGER NCOADDS_ARRAY      !number of coadds (needed to compute errors)
        INTEGER NBOX_CENTROID                         !box size to fit centroid
        INTEGER NSIMUL,NSIMUL_                           !number of simulations
        INTEGER IMODEG2D                              !fit mode for 2D gaussian
        INTEGER NBUFF_R,NBUFF_G,NBUFF_B            !buffers for RGB *.ppm files
        INTEGER NFITSFILES             !number of fits files to measure offsets
        INTEGER NCOLASC1(NMAXBUFF),NCOLASC2(NMAXBUFF)    !plot marks from ASCII
        INTEGER ASCCOLOR(NMAXBUFF),ASCLWIDTH(NMAXBUFF)   !plot marks from ASCII
        INTEGER ASCSYMB(NMAXBUFF)                        !plot marks from ASCII
        INTEGER ASCNUMBER(NMAXBUFF)              !plot symbol number from ASCII
        INTEGER NIARG                                !number of input arguments
        REAL XC,YC,XC_,YC_,XC__,YC__                            !mouse location
        REAL XOFFSET(9),YOFFSET(9)                            !measured offsets
        REAL BG,FG                              !image bacground and foreground
        REAL BG_,FG_
        REAL Z1,Z2                                       !ZSCALE al estilo Iraf
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)                        !current image
!delete REAL IMAGEN_(NXMAX,NYMAX)
        REAL FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX  !basic statistic of current region
        REAL FMEAN_,FSIGMA_
        REAL XCUTX(NXYMAX),YCUTX(NXYMAX)     !X cut
        REAL XP(NXYMAX),YP(NXYMAX)               !generic matrix for X or Y cut 
        REAL XSTACK(NXYMAX),YSTACK(NXYMAX)      !coordinates of points in stack
        REAL COEFFBL00(20),COEFFBA00(20)
        REAL FSHIFT                                        !Y shift of BL lines
        REAL COEFF(20)                                 !polynomial coefficients
        REAL CWV(20)         !polynomial coefficients of wavelength calibration
        REAL XMINFIT,XMAXFIT         !limits to perform fit of refined BL lines
        REAL AIJ(NECUAMAX),BIJ(NECUAMAX)
        REAL AIJ_(NECUAMAX),BIJ_(NECUAMAX)
        REAL U,V
        REAL TSIGMA                      !times sigma to exclude points in fits
        REAL GAIN_ARRAY                          !array gain (in electrons/ADU)
        REAL RNOISE_ARRAY                         !array readout noise (in ADU)
        REAL YRMSTOL                                              !for DOWNHILL
        REAL XMIN,XMAX,YMIN,YMAX
        REAL FRGB_R,FRGB_G,FRGB_B   !factors to scale images in RGB composition
        REAL FLOGRGB                              !factor for logarithmic scale
        REAL X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE     !coefficients of centroid fit
        REAL EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE !errors
        REAL XDUMMY,YDUMMY                               !dummy float variables
        REAL PHOTFLAM(NMAXBUFF),PHOTZPT(NMAXBUFF),EXPTIME(NMAXBUFF)   !HST flux
        REAL XOFFSET_ORIGEN,YOFFSET_ORIGEN           !origen para medir offsets
        REAL STWV,DISP,AIRMASS,TIMEXPOS
        REAL MASKVALUE                        !value to be masked in statistics
        REAL CRPIX1(NMAXBUFF),CRVAL1(NMAXBUFF),CDELT1(NMAXBUFF)  !wavel. calib.
        REAL ASCCHEIGHT(NMAXBUFF)                        !plot marks from ASCII
        REAL XX0(1),YY0(1)
        CHARACTER*1 CH                            !mouse button or keyboard key
        CHARACTER*1 CDUM                             !dum character*1 variables
        CHARACTER*1 CSAVE               !choose between FITS or REDUCEME format
        CHARACTER*1 CINERR             !determina el tipo de fichero de errores
        CHARACTER*1 CSAVEERR       !determina si salvamos el fichero de errores
        CHARACTER*1 CSURE                      !confirmation for end of program
        CHARACTER*1 CREFINE                      !refine BA or BL line location
        CHARACTER*1 CUPDATE                              !update BA or BL lines
        CHARACTER*1 CRGBSATUR      !correction of saturation in RGB composition
        CHARACTER*1 CENHANCE_RB        !enhance red and blue in RGB composition
        CHARACTER*1 CLOGRGB        !use of logarithmic scale in RGB composition
        CHARACTER*1 COFFSET !indicates whether we measure offsets with centroid
        CHARACTER*1 CSPECIAL               !option in menu of special utilities
        CHARACTER*1 COFFSETS       !mode to compute offsets with list of images
        CHARACTER*20 CIXC,CWXC,CIYC,CSIGNAL!character strings to display cursor
        CHARACTER*50 CDUMMY                             !dum character variable
        CHARACTER*255 TTER                                    !graphics display
        CHARACTER*250 CLINBUT               !character string with accelerators
        CHARACTER*255 INFILE_,INFILEBOX9,ERRFILE       !generic input file name
        CHARACTER*255 OFFFILE                          !generic input file name
        CHARACTER*255 INFILE(NMAXBUFF)                         !input file name
        CHARACTER*255 WAVEFILE  !file name of polynomial wavelength calibration
        CHARACTER*255 OUTFILE                         !generic output file name
        CHARACTER*255 HDRFILE                         !generic output file name
        CHARACTER*255 NEWDEV      !generic output file name for postscript plot
        CHARACTER*255 FILELISTIN             !file name with list of FITS files
        CHARACTER*255 FILELISTOUT            !file name with list of FITS files
        CHARACTER*255 FILEOFFSET
        CHARACTER*255 DS9REGFILE(NMAXBUFF)                    !ds9 region files
        CHARACTER*255 ASCREGFILE(NMAXBUFF)                    !ds9 region files
        CHARACTER*255 CARG(NMAXBUFF/2)
        LOGICAL LEXIT                              !govern the main button loop
        LOGICAL LBEXIST                !button selected with keyboard is active
        LOGICAL LOGFILE,LOGFILERR             !govern the reading of a new file
        LOGICAL LNULL(NXMAX,NYMAX,NMAXBUFF),ANYNULL !NaN, Infty, etc. in IMAGEN
        LOGICAL LFIRSTPLOT                                  !next is first plot
        LOGICAL LSTACK                                          !activate stack
        LOGICAL LASK                  !generic logical variable to govern loops
        LOGICAL LOK                      !if .TRUE. a routine has done properly
        LOGICAL LBOX9    !if .TRUE. INFILE_ is name of list with frames of BOX9
        LOGICAL LECHO         !if .TRUE., the I/O functions echo the input data
        LOGICAL LMIDE                       !if .TRUE. measure offsets in box-9
        LOGICAL LREPEAT
        LOGICAL L_PHOTFLAM(NMAXBUFF),L_PHOTZPT(NMAXBUFF)              !HST flux
        LOGICAL L_EXPTIME(NMAXBUFF)                                   !HST flux
        LOGICAL LOVERCUTS
        LOGICAL LDS9REG(NMAXBUFF)               !if .TRUE. overplot ds9 regions
        LOGICAL LASCREG(NMAXBUFF)      !if .TRUE. overplot X,Y from binary FITS
        LOGICAL LWAVECAL(NMAXBUFF) !if .TRUE., image has wavelength calibration
! common blocks
!delete COMMON/BLKIMAGEN1/IMAGEN             !imagen FITS leida en formato REAL
!delete COMMON/BLKIMAGEN1_/IMAGEN_              !es global para ahorrar memoria
        COMMON/BLKIMAGEN2/NCBUFF                      !numero del buffer actual
        COMMON/BLKLNULL/LNULL,ANYNULL   !mascara que indica si existen NaN, etc
        COMMON/BLKNAXIS/NAXIS                                      !dimensiones
        COMMON/BLKNAXISFRAME/NAXISFRAME
        COMMON/BLKNSIZEFB9/NSIZEFB9
        COMMON/BLKNFRAMES/NFRAMES                  !real number of frames/image
        COMMON/BLKJUST/JUST
        COMMON/BLKBGFG/BG,FG
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKESTADISTICA/FMEAN,FSIGMA,FMEDIAN,FMIN,FMAX
        COMMON/BLKINFILE/INFILE
        COMMON/BLKBOUND1/NLINBL,NDEGBL,NDEGBL00
        COMMON/BLKBOUND1B/NLINBA,NDEGBA,NDEGBA00
        COMMON/BLKBOUND2/COEFFBL
        COMMON/BLKBOUND2B/COEFFBA
        COMMON/BLKBOUND3/LBOUNDL,LBOUNDA
        COMMON/BLKBOUND4/XMINBL,YMINBL,XMAXBL,YMAXBL
        COMMON/BLKBOUND5/XMINBA,YMINBA,XMAXBA,YMAXBA
        COMMON/BLKMAPPING0/NMAP
        COMMON/BLKMAPPING1/AIJ,BIJ
        COMMON/BLKMAPPING2/AIJ_,BIJ_
        COMMON/BLKMAPPING3/SXMINGRID,SYMINGRID,SXMAXGRID,SYMAXGRID
        COMMON/BLKMAPPING4/SXMINEXTG,SYMINEXTG,SXMAXEXTG,SYMAXEXTG
        COMMON/BLKDEFAULTS1/NYFINDAL,NINLINE
        COMMON/BLKDEFAULTS2/NWIDTH
        COMMON/BLKDEFAULTS3/TSIGMA
        COMMON/BLKDEFAULTS4/NGRID
        COMMON/BLKDEFAULTS5/IRESAMPL
        COMMON/BLKDEFAULTS6/GAIN_ARRAY
        COMMON/BLKDEFAULTS7/RNOISE_ARRAY
        COMMON/BLKDEFAULTS8/NBOX_CENTROID
        COMMON/BLKDEFAULTS9/YRMSTOL
        COMMON/BLKDEFAULTS10/NSIMUL
        COMMON/BLKDEFAULTS11/IMODEG2D
        COMMON/BLKDEFAULTS12/XOFFSET_ORIGEN,YOFFSET_ORIGEN
        COMMON/BLKRGBBUFF/NBUFF_R,NBUFF_G,NBUFF_B
        COMMON/BLKRGBFACT/FRGB_R,FRGB_G,FRGB_B,FLOGRGB
        COMMON/BLKRGBSATUR/CRGBSATUR,CENHANCE_RB,CLOGRGB
        COMMON/BLKLECHO/LECHO
        COMMON/BLK_HSTFLUX1/L_PHOTFLAM,L_PHOTZPT,L_EXPTIME !keywords with info
        COMMON/BLK_HSTFLUX2/PHOTFLAM,PHOTZPT,EXPTIME !concerning HST flux cal.
        COMMON/BLKREDUCEME/STWV,DISP,AIRMASS,TIMEXPOS
        COMMON/BLKLOVERCUTS/LOVERCUTS
        COMMON/BLKLDS9REG/LDS9REG
        COMMON/BLKLASCREG/LASCREG
        COMMON/BLKDS9REGFILE/DS9REGFILE
        COMMON/BLKASCREGFILE/ASCREGFILE
        COMMON/BLKASCCOLS1/NCOLASC1,NCOLASC2,ASCCOLOR,ASCLWIDTH,ASCSYMB
        COMMON/BLKASCCOLS2/ASCCHEIGHT,ASCNUMBER
        COMMON/BLKWAVECAL1/CRPIX1,CRVAL1,CDELT1
        COMMON/BLKWAVECAL2/LWAVECAL
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        CALL Initialize_Dynamic_Array_IMAGEN
        CALL Initialize_Dynamic_Array_IMAGEN_
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Note: the pattern of the frames in box-9 is the following:
!       6 9 4
!       3 1 7
!       8 5 2
! offsets of each frame in the 768x768 composite mosaic
        DATA (DI(K),DJ(K),K=1,NBOXMAX) / &
!!!     +   256,256,  !1
!!!     +   000,512,  !2
!!!     +   256,000,  !3
!!!     +   512,512,  !4
!!!     +   000,256,  !5
!!!     +   512,000,  !6
!!!     +   256,512,  !7
!!!     +   000,000,  !8
!!!     +   512,256/  !9
         1,1, & !1
         0,2, & !2
         1,0, & !3
         2,2, & !4
         0,1, & !5
         2,0, & !6
         1,2, & !7
         0,0, & !8
         2,1/  !9
        LBOX9=.FALSE.
!
        I1=0 !evita WARNING de compilacion
        J1=0 !evita WARNING de compilacion
        NNX1_=0 !evita WARNING de compilacion
        NNY1_=0 !evita WARNING de compilacion
! Defaults for symbols at (X,Y) coordinates given in external ASCII file
        DO I=1,NMAXBUFF
          NCOLASC1(I)=1
          NCOLASC2(I)=2
          ASCSYMB(I)=25
          ASCCOLOR(I)=3
          ASCLWIDTH(I)=1
          ASCCHEIGHT(I)=1.0
          ASCNUMBER(I)=0   !0=no, 1=yes
          ASCREGFILE(I)='none'
        END DO
!------------------------------------------------------------------------------
        WRITE(*,101) 'Welcome to NIRSPEC'
        WRITE(*,101) 'This is version 5.0'
!
        NIARG=IARGC()
        IF(NIARG.GT.NMAXBUFF/2)THEN
          WRITE(*,100) 'Number of input arguments: '
          WRITE(*,*) NIARG
          WRITE(*,100) 'Maximum number allowed...: '
          WRITE(*,*) NMAXBUFF/2
          WRITE(*,101) 'FATAL ERROR: number of input arguments is too large'
          CALL Deallocate_Array_IMAGEN
          CALL Deallocate_Array_IMAGEN_
          STOP
        END IF
!       DO I=1,NIARG
!         CALL GETARG(I, CARG(I))
!         L1=TRUELEN(CARG(I))
!         WRITE(*,*) '-->'//CARG(I)(1:L1)//'<--'
!       END DO
!       STOP
!------------------------------------------------------------------------------
        NXYMAX_ = AMAX0(NXMAX, NYMAX)
        IF(NXYMAX.NE.NXYMAX_)THEN
          WRITE(*,100) 'NXMAX...: '
          WRITE(*,*) NXMAX
          WRITE(*,100) 'NYMAX...: '
          WRITE(*,*) NYMAX
          WRITE(*,100) 'NXYMAX..: '
          WRITE(*,*) NXYMAX
          WRITE(*,101) 'FATAL ERROR: NXYMAX must be set to the MAX(NXMAX,NYMAX) value in configure.ac'
          CALL Deallocate_Array_IMAGEN
          CALL Deallocate_Array_IMAGEN_
          STOP
        END IF
!------------------------------------------------------------------------------
        NCBUFF=1
        DO K=1,NMAXBUFF
          NAXIS(1,K)=0
          NAXIS(2,K)=0
          NFRAMES(K)=0
          INFILE(K)='[undefined]'
          L_PHOTFLAM(K)=.FALSE.
          L_PHOTZPT(K)=.FALSE.
          L_EXPTIME(K)=.FALSE.
          PHOTFLAM(K)=0.0
          PHOTZPT(K)=0.0
          EXPTIME(K)=0.0
          LDS9REG(K)=.FALSE.
          LASCREG(K)=.FALSE.
        END DO
!
        LFIRSTPLOT=.TRUE.
        LBOUNDARY=.FALSE.
        LBOUNDL=.FALSE.
        LBOUNDA=.FALSE.
        LPBOUND=.FALSE.
        LMAPPING=.FALSE.
        LPMAPP=.FALSE.
        LPWAVE=.FALSE.
        LOVERCUTS=.FALSE.
        JUST=0
!
        MODECUT=1 !single
        LSTACK=.FALSE.
        NSTACK=0
!
        NCOLOR1=-3-1 !load, Save,...
        NCOLOR2=-2-1 !zoom, Whole,...
        NCOLOR3=-4-1 !histogram...
        NCOLOR4=-5-1 !Plot commands...
        NCOLOR5=-7-1 !detect sky lines...
        NCOLOR6=-6-1 !image operations...
!
        CINERR='0'
        CSAVE='f'
        CSAVEERR='n'
        NBOX=0
!!!     IROTATE=-90 !default for NIRSPEC FITS images
        IROTATE=0
!
        NYFINDAL=40     !Y grid looking for spectral lines
        NINLINE=25      !no. of peaks/line found (max. NYFINDAL)
        NWIDTH=11       !no. of pixels to search for peaks (odd, min. 3)
        TSIGMA=3.0      !times sigma to exclude points in fits
        NGRID=20        !grid width (pixels) to plot mapping
        IRESAMPL=1      !resampling method (1=approximate,2=accurate)
        GAIN_ARRAY=5.   !array gain (in electrons/ADU)
        RNOISE_ARRAY=5. !array readout noise (in ADU)
        NCOADDS_ARRAY=1 !number of coadds
        NBOX_CENTROID=15!box size to fit centroid
        YRMSTOL=1.E-5   !for DOWNHILL
        NSIMUL=0        !number of simulations
        IMODEG2D=0      !fit mode for 2D gaussian
        NBUFF_R=-1      !undefined
        NBUFF_G=-1      !undefined
        NBUFF_B=-1      !undefined
        FRGB_R=1.0
        FRGB_G=1.0
        FRGB_B=1.0
        CRGBSATUR='n'
        CENHANCE_RB='n'
        CLOGRGB='n'
        FLOGRGB=2.0
        NBUFFBOX9=0
        XOFFSET_ORIGEN=0.0
        YOFFSET_ORIGEN=0.0
        COFFSETS='u'
! Por si usamos formato REDUCEME
        STWV=0.
        DISP=0.
        AIRMASS=0.
        TIMEXPOS=0.
!
        MASKVALUE=0.0 !valor a ignorar en la estadistica (en caso necesario)
!------------------------------------------------------------------------------
! Open graphic display and establish button and plot settings
        TTER='/XSERVE'
        TTER=READC('Graphic device',TTER,'@')
        IF(TTER(1:1).EQ.'@')THEN
          TTER=TTER(2:)
          LECHO=.TRUE.
          CALL RPGBEGOK(TTER,1)
        ELSE
          LECHO=.FALSE.
          CALL RPGBEGOK(TTER,0)
        END IF
        CALL PALETTE(3)
        CALL PGQID(ID)
        CALL BUTTSPR(0.440,0.965,0.050,0.750)
        CALL BUTTSBR(0.000,1.000,0.000,1.000)
        CALL BUTTSXB(10)
        CALL BUTTSYB(25)
        CALL BUTTSCF(1)
        CALL BUTTSCH(0.8)
!------------------------------------------------------------------------------
        CALL BUTTON(1,'[l]oad',0)
        CALL BUTTON(1,'[l]oad',NCOLOR1)
        CALL BUTTON(2,'[s]ave',0)
        CALL BUTTON(2,'[s]ave',3)
        CALL BUTTON(11,'[q]uit',0)
        CALL BUTTON(11,'[q]uit',NCOLOR1)
        CALL BUTTON(21,'[?]buffers',0)
        CALL BUTTON(22,'postscript',0)
        CALL BUTTON(22,'postscript',3)
        CALL BUTTON(23,'RGB *.ppm',0)
        CALL BUTTON(23,'RGB *.ppm',3)
        CALL BUTTON(12,'[i]math',0)
        CALL BUTTON(12,'[i]math',3)
!
        CALL BUTTON(3,'[b]oundary',0)
        CALL BUTTON(3,'[b]oundary',3)
        CALL BUTTON(13,'plot bou.',0)
        CALL BUTTON(13,'plot bou.',3)
!
        CALL BUTTON(5,'[z]oom',0)
        CALL BUTTON(5,'[z]oom',3)
        CALL BUTTON(6,'[w]hole',0)
        CALL BUTTON(6,'[w]hole',3)
        IF(JUST.EQ.1)THEN
          CALL BUTTON(7,'x[=]y',0)
          CALL BUTTON(7,'x[=]y',3)
        ELSE
          CALL BUTTON(7,'  x![=]y',0)
          CALL BUTTON(7,'  x![=]y',3)
        END IF
        IF(MODECUT.EQ.1)THEN
          CALL BUTTON(15,'single',0)
          CALL BUTTON(15,'single',3)
        ELSEIF(MODECUT.EQ.2)THEN
          CALL BUTTON(15,'region',0)
          CALL BUTTON(15,'region',3)
        ELSE
          CALL BUTTON(15,'all',0)
          CALL BUTTON(15,'all',3)
        END IF
        CALL BUTTON(16,'[x] cut',0)
        CALL BUTTON(16,'[x] cut',3)
        CALL BUTTON(17,'[y] cut',0)
        CALL BUTTON(17,'[y] cut',3)
!
        CALL BUTTON(4,'find lines 1',0)
        CALL BUTTON(4,'find lines 1',3)
        CALL BUTTON(14,'find lines 2',0)
        CALL BUTTON(14,'find lines 2',3)
        CALL BUTTON(8,'[m]apping',0)
        CALL BUTTON(8,'[m]apping',3)
        CALL BUTTON(8,'[m]apping',0)
        CALL BUTTON(8,'[m]apping',3)
        CALL BUTTON(9,'wa[v]el.cal.',0)
        CALL BUTTON(9,'wa[v]el.cal.',3)
        CALL BUTTON(18,'[p]lot mapp.',0)
        CALL BUTTON(18,'[p]lot mapp.',3)
        CALL BUTTON(19,'plot wave.',0)
        CALL BUTTON(19,'plot wave.',3)
!
        CALL BUTTON(10,'stack 0',0)
        IF(.NOT.LSTACK)THEN
          CALL BUTTON(10,'stack 0',3)
        END IF
!
        CALL BUTTON(20,'s[Tt]atistic',0)
        CALL BUTTON(20,'s[Tt]atistic',3)
!
        CALL BUTTON(24,'[Cc]entroid',0)
        CALL BUTTON(24,'[Cc]entroid',3)
!
        CALL BUTTON(25,'[r]ectify',0)
        CALL BUTTON(25,'[r]ectify',3)
        CALL BUTTON(26,'s[k]y sub.',0)
        CALL BUTTON(26,'s[k]y sub.',3)
        CALL BUTTON(27,'X-[e]xtract',0)
        CALL BUTTON(27,'X-[e]xtract',3)
        CALL BUTTON(28,'[u]nrectify',0)
        CALL BUTTON(28,'[u]nrectify',3)
        CALL BUTTON(29,'[.]special',0)
        CALL BUTTON(29,'[.]special',3)
        CALL BUTTON(30,'offsets',0)
!
        CALL BUTTON(31,'resize ima',0)
        CALL BUTTON(31,'resize ima',3)
        CALL BUTTON(32,'combine',0)
        CALL BUTTON(32,'combine',3)
        CALL BUTTON(33,'shiftb9',0)
        CALL BUTTON(33,'shiftb9',3)
!
        CALL BUTTON(34,'box[9]oper',0)
        CALL BUTTON(34,'box[9]oper',3)
        CALL BUTTON(44,'[o]ffsets',0)
        CALL BUTTON(44,'[o]ffsets',3)
!
        CALL BUTTON(41,'zoom',0)
        CALL BUTTON(41,'zoom',3)
        CALL BUTTON(42,'whole',0)
        CALL BUTTON(42,'whole',3)
        CALL BUTTON(43,'min,max',0)
        CALL BUTTON(43,'min,max',3)
!
        CALL BUTTON(35,'#[1]',0)
        CALL BUTTON(36,'#[2]',0)
        CALL BUTTON(37,'#[3]',0)
        CALL BUTTON(38,'#[4]',0)
        CALL BUTTON(39,'#[5]',0)
        CALL BUTTON(40,'#[6]',0)
        CALL BUTTON(45,'err #1',0)
        CALL BUTTON(46,'err #2',0)
        CALL BUTTON(47,'err #3',0)
        CALL BUTTON(48,'err #4',0)
        CALL BUTTON(49,'err #5',0)
        CALL BUTTON(50,'err #6',0)
!
        CALL BUTTON(161,'zoom',0)
        CALL BUTTON(161,'zoom',3)
        CALL BUTTON(162,'min[,]max',0)
        CALL BUTTON(162,'min[,]max',3)
        CALL BUTTON(163,'z1[/]z2',0)
        CALL BUTTON(163,'z1[/]z2',3)
        CALL BUTTON(164,'BG[:]FG',0)
        CALL BUTTON(164,'BG[:]FG',3)
!------------------------------------------------------------------------------
! Load images indicated as arguments in the command line
        IF(NIARG.GT.0)THEN
          DO NEWBUFF=1,NIARG
            CALL GETARG(NEWBUFF, CARG(NEWBUFF))
            L1=TRUELEN(CARG(NEWBUFF))
            INFILE_=CARG(NEWBUFF)(1:L1)
            INQUIRE(FILE=INFILE_,EXIST=LOGFILE)
            IF(LOGFILE)THEN
              WRITE(*,100) 'Reading file: '
              WRITE(*,101) INFILE_
              CALL SLEEFITS(INFILE_,.FALSE.,IROTATE,NEWBUFF,LBOX9,NEWBUFF)
              WRITE(*,101) ' ...OK!'
              IF(ANYNULL) THEN
                WRITE(*,101) '***WARNING***'
                WRITE(*,101) '=> ANYNULL=.TRUE.!'
                !CALL PGEND
                !STOP
              END IF
              NFRAMES(NEWBUFF)=1
              NSIZEFB9(NEWBUFF)=0
              NAXISFRAME(1,1,NEWBUFF)=0
              NAXISFRAME(2,1,NEWBUFF)=0
            ELSE
              WRITE(*,101) 'FATAL ERROR: the following file does not exist:'
              WRITE(*,101) INFILE_
              CALL Deallocate_Array_IMAGEN
              CALL Deallocate_Array_IMAGEN_
              STOP
            END IF
            NCBUFF=NEWBUFF
            INFILE(NCBUFF)=INFILE_
            DO K=1,NMAXBUFF
              IF(K.LE.NMAXBUFF/2)THEN
                WRITE(CDUMMY,'(A2,I1,A1)') '#[',K,']'
                CALL RMBLANK(CDUMMY,CDUMMY,L)
                IF(K.EQ.NCBUFF)THEN
                  CALL BUTTON(34+K,CDUMMY(1:L),1)
                ELSE
                  CALL BUTTON(34+K,CDUMMY(1:L),0)
                END IF
              ELSE
                WRITE(CDUMMY,'(A5,I1)') 'err #',K-NMAXBUFF/2
                L=TRUELEN(CDUMMY)
                IF(K.EQ.NCBUFF)THEN
                  CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),1)
                ELSE
                  CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),0)
                END IF
              END IF
            END DO
            NX1=1
            NX2=NAXIS(1,NCBUFF)
            NY1=1
            NY2=NAXIS(2,NCBUFF)
            CALL STATISTIC(NCBUFF,NX1,NX2,NY1,NY2,.FALSE.,.TRUE.,.TRUE.,0.0,.FALSE.)
            IF(LFIRSTPLOT)THEN
              CALL ACTIVEBUT(JUST,MODECUT)
              IF(FSIGMA.GT.0.0)THEN
                BG=FMEAN-5.*FSIGMA
                FG=FMEAN+5.*FSIGMA
              ELSE
                BG=FMEAN-1.0
                FG=FMEAN+1.0
              END IF
              CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
            ELSE
              CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
            END IF
            IF(LFIRSTPLOT)THEN
              CALL HISTOGRAM(NCBUFF)
              LFIRSTPLOT=.FALSE.
            END IF
            NY_CUT=NY2/2
            IF(NY_CUT.EQ.0) NY_CUT=1
            DO J=1,NAXIS(1,NCBUFF)
              XCUTX(J)=REAL(J)
              YCUTX(J)=IMAGEN(J,NY_CUT,NCBUFF)
            END DO
            WRITE(CDUMMY,'(A1,I10,A1,I10,A1)') '[',NY_CUT,',',NY_CUT,']'
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL SUBPLOT(NAXIS(1,NCBUFF),1,NAXIS(1,NCBUFF),XCUTX,YCUTX,XCUTX,YCUTX,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
             'x axis','signal',CDUMMY(1:L),NCBUFF,201,1.0)
            CALL BUTTON(51,'over=FALSE',0)
            CALL BUTTON(51,'over=FALSE',NCOLOR4)
          END DO
        END IF
!------------------------------------------------------------------------------
! ***MAIN LOOP***
!------------------------------------------------------------------------------
        LEXIT=.FALSE.
        DO WHILE(.NOT.LEXIT)
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
! change uppercase by lowercase when appropriate
          IF(CH.EQ.'Z') CH='z'
          IF(CH.EQ.'W') CH='w'
! look for accelerators
          IF(ICHAR(CH).EQ.9)THEN !tabulador
            NBLOCAL=15
          ELSEIF(CH.EQ.'Q')THEN
            NBLOCAL=11
          ELSEIF((CH.EQ.'X').AND.LSTACK)THEN
            NBLOCAL=10
          ELSEIF(CH.EQ.'$')THEN !cambia entre imagen de errores y datos
            IF(NCBUFF.LE.NMAXBUFF/2)THEN
              NBLOCAL=44+NCBUFF
            ELSE
              NBLOCAL=28+NCBUFF
            END IF
          ELSEIF(CH.EQ.'C')THEN !centering
            NBLOCAL=24
          ELSEIF(CH.EQ.'T')THEN !estadistica con el ultimo taman~o de caja
            NBLOCAL=20
          ELSE
!.....................0000000001
!.....................1234567890
            CLINBUT=('lsb zw=mv '// & !00
                     'qi   xyp t'// & !01
                     '?  crkeu. '// & !02
                     '   9123456'// & !03
                     '   o      '// & !04
                     '          '// & !05
                     '          '// & !06
                     '          '// & !07
                     '          '// & !08
                     '          '// & !09
                     '          '// & !10
                     '          '// & !11
                     '          '// & !12
                     '          '// & !13
                     '          '// & !14
                     '          '// & !15
                     ' ,/:      '// & !16
                     '          '// & !17
                     '          '// & !18
                     '          '// & !19
                     '          '// & !20
                     '          '// & !21
                     '          '// & !22
                     '          '// & !23
                     '          ')  !24
            NBLOCAL=INDEX(CLINBUT,CH)
          END IF
          IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
            CALL BUTTQEX(NBLOCAL,LBEXIST)
            IF(LBEXIST) NB=NBLOCAL
          END IF
!------------------------------------------------------------------------------
!!!     print*,'nb=',nb
          IF(NB.EQ.0)THEN
!..............................................................................
            IF(CH.EQ.'A')THEN
              IF(.NOT.LFIRSTPLOT)THEN
                IXC=NINT(XC)
                IYC=NINT(YC)
                IF((IXC.GE.NX1).AND.(IXC.LE.NX2).AND.(IYC.GE.NY1).AND.(IYC.LE.NY2))THEN
                  WRITE(CIXC,*) XC
                  CALL RMBLANK(CIXC,CIXC,L1)
                  IF(LWAVECAL(NCBUFF))THEN
                    WRITE(CWXC,*) (XC-CRPIX1(NCBUFF))*CDELT1(NCBUFF)+CRVAL1(NCBUFF)
                    CALL RMBLANK(CWXC,CWXC,L2)
                  END IF
                  WRITE(CIYC,*) YC
                  CALL RMBLANK(CIYC,CIYC,L3)
                  WRITE(CSIGNAL,*) IMAGEN(IXC,IYC,NCBUFF)
                  CALL RMBLANK(CSIGNAL,CSIGNAL,L4)
                  IF(LWAVECAL(NCBUFF))THEN
                    WRITE(*,101) '=> Cursor at X='//CIXC(1:L1)//', WV='//CWXC(1:L2)//', Y='//CIYC(1:L3)//', signal='//CSIGNAL(1:L4)
                  ELSE
                    WRITE(*,101) '=> Cursor at X='//CIXC(1:L1)//', Y='//CIYC(1:L3)//', signal='//CSIGNAL(1:L4)
                  END IF
                  IF(LMAPPING)THEN
                    CALL FMAP(NMAP,AIJ,BIJ,XC,YC,U,V)
                    WRITE(CIXC,*) U
                    CALL RMBLANK(CIXC,CIXC,L1)
                    WRITE(CIYC,*) V
                    CALL RMBLANK(CIYC,CIYC,L2)
                    WRITE(*,101) '=> Sx='//CIXC(1:L1)//', Sy='//CIYC(1:L2)
                    IF(LPWAVE)THEN
                      WRITE(*,100) '=> Wavelength: '
                      WRITE(*,*) FPOLY(NDEGWV,CWV,U+0.5+REAL(SXMINEXTG))
                    END IF
                  END IF
                  IF(LSTACK)THEN
                    IF(NSTACK.EQ.NXYMAX)THEN
                      WRITE(*,101) '***ERROR***'
                      WRITE(*,100) '=> stack is full'
                      WRITE(*,100) ' (press <CR> to continue...)'
                      READ(*,*)
                    ELSE
                      CALL PGSCI(2)
                      !usamos un array unidimensional porque el compilador
                      !gfortran-mp-10 da error al usar un escalar en lugar
                      !de una matriz
                      XX0(1)=REAL(IXC)
                      YY0(1)=REAL(IYC)
                      CALL PGPOINT(1,XX0,YY0,21)
                      CALL PGSCI(1)
                      WRITE(*,101) 'Please, confirm this point...'
                      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
                      IF(CH.EQ.'A')THEN
                        WRITE(*,101) 'OK! Point added to stack.'
                        WRITE(*,*)
                        NSTACK=NSTACK+1
                        XSTACK(NSTACK)=REAL(IXC)
                        YSTACK(NSTACK)=REAL(IYC)
                        CALL PGSCI(3)
                        CALL PGPOINT(1,XSTACK(NSTACK),YSTACK(NSTACK),21)
                        CALL PGSCI(1)
                        WRITE(CDUMMY,*) NSTACK
                        CALL RMBLANK(CDUMMY,CDUMMY,L1)
                        CALL BUTTON(10,'stack '//CDUMMY(1:L1),0)
                        CALL BUTTON(10,'stack '//CDUMMY(1:L1),1)
                      ELSE
                        WRITE(*,101) 'WARNING: point canceled.'
                        WRITE(*,*)
                      END IF
                    END IF
                  END IF
                END IF
              END IF
!..............................................................................
            ELSE
              print*,'ichar: ',ichar(ch)
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.1)THEN
            CALL BUTTON(NB,'[l]oad',5)
            LOGFILE=.FALSE.
            LOGFILERR=.FALSE.
            DO WHILE(.NOT.LOGFILE)
              INFILE_=READC('Input FITS file name (none=EXIT)','*.*fit* FITS*/*.fits *.b9','@')
              IF((INDEX(INFILE_,'*').NE.0).OR.(INDEX(INFILE_,'?').NE.0))THEN
                L1=TRUEBEG(INFILE_)
                L2=TRUELEN(INFILE_)
                ISYSTEM=SYSTEMFUNCTION('ls '//INFILE_(L1:L2))
              ELSEIF(INFILE_.EQ.'none')THEN
                LOGFILE=.TRUE.
                LOGFILERR=.TRUE.
              ELSE
                L2=TRUELEN(INFILE_)
                IF(INFILE_(L2-2:L2).EQ.'.b9')THEN
                  LBOX9=.TRUE.
                ELSE
                  LBOX9=.FALSE.
                END IF
                INQUIRE(FILE=INFILE_,EXIST=LOGFILE)
                IF(.NOT.LOGFILE)THEN
                  L1=TRUEBEG(INFILE_)
                  L2=TRUELEN(INFILE_)
                  INQUIRE(FILE=INFILE_(L1:L2)//'.fits',EXIST=LOGFILE)
                  IF(LOGFILE)THEN
                    INFILE_=INFILE_(L1:L2)//'.fits'
                    WRITE(*,100) '>>> Assuming fits extension: '
                    WRITE(*,101) INFILE_(L1:L2)//'.fits'
                  ELSE
                    CALL LEEONEFILE(INFILE_,LOK)
                    IF(LOK)THEN
                      L1=TRUEBEG(INFILE_)
                      L2=TRUELEN(INFILE_)
                      WRITE(*,100) '>>> Using unambiguous file name: '
                      WRITE(*,101) INFILE_(L1:L2)
                      LOGFILE=.TRUE.
                    ELSE
                      WRITE(*,101) '***ERROR***'
                      WRITE(*,100) '=> This file does not exist.'
                      WRITE(*,100) ' Try again. (press <CR>...)'
                      READ(*,*)
                    END IF
                  END IF
                END IF
              END IF
            END DO
!..............................................................................
            IF(.NOT.LOGFILERR)THEN                 !the file name is not 'none'
              LASK=.TRUE.
              DO WHILE(LASK)
                IF(LBOX9)THEN
                  IROTATE=0
                ELSE
                  WRITE(CDUMMY,*) IROTATE
                  IROTATE=READI('Rotation angle to be applied',CDUMMY)
                END IF
                LASK=(IROTATE.NE.0).AND.(IABS(IROTATE).NE.90).AND.(IABS(IROTATE).NE.180).AND.(IABS(IROTATE).NE.270)
                IF(LASK)THEN
                  WRITE(*,101) '***ERROR***'
                  WRITE(*,100) '=> Invalid angle.'
                  WRITE(*,101) ' Only 0, +/-90, +/-180 +/-270'
                  WRITE(*,101) 'deg. are valid options.'
                  WRITE(*,100) '(press <CR> to continue...)'
                  READ(*,*)
                END IF
              END DO
              IF(LFIRSTPLOT)THEN
                NEWBUFF=1
              ELSE
                NEWBUFF=NCBUFF+1
                IF(NEWBUFF.GT.NMAXBUFF/2) NEWBUFF=1
              END IF
              WRITE(CDUMMY,*) NEWBUFF
              NEWBUFF=READILIM('New buffer #',CDUMMY,1,NMAXBUFF)
              IF(NEWBUFF.LE.NMAXBUFF/2)THEN
                WRITE(*,101) '* ERROR frame:'
                WRITE(*,101) '(0) initialize to zero'
                WRITE(*,101) '(g) compute from GAIN and RN'
                WRITE(*,101) '(f) additional file'
                IF(CINERR.EQ.'x') CINERR='0'
                CINERR(1:1)=READC('Option',CINERR,'0gf')
                IF(CINERR.EQ.'g')THEN
                  WRITE(CDUMMY,*) GAIN_ARRAY
                  GAIN_ARRAY=READF(  'Gain........(e/ADU)',CDUMMY)
                  WRITE(CDUMMY,*) RNOISE_ARRAY
                  RNOISE_ARRAY=READF('Readout noise (ADU)',CDUMMY)
                  LASK=.TRUE.
                  DO WHILE(LASK)
                    WRITE(CDUMMY,*) NCOADDS_ARRAY
                    NCOADDS_ARRAY=READI('Number of coadds...',CDUMMY)
                    LASK=(NCOADDS_ARRAY.LE.0)
                    IF(LASK)THEN
                      WRITE(*,101) '***ERROR***'
                      WRITE(*,101) '=> Number of coadds must be > 0'
                      WRITE(*,100) '(press <CR> to continue...)'
                      READ(*,*)
                    END IF
                  END DO
                END IF
              ELSE
                CINERR='x'
              END IF
              IF(LBOX9)THEN
                CALL MAXDIMB9(NEWBUFF,INFILE_,NSIZEFB9(NEWBUFF))
                WRITE(*,100) '> Maximum size of frame in box-9: '
                WRITE(*,*) NSIZEFB9(NEWBUFF)
                NAXIS(1,NEWBUFF)=3*NSIZEFB9(NEWBUFF)
                NAXIS(2,NEWBUFF)=3*NSIZEFB9(NEWBUFF)
                IF(NAXIS(1,NEWBUFF).GT.NXMAX)THEN
                  WRITE(*,100) 'NAXIS(1), NXMAX: '
                  WRITE(*,*) NAXIS(1,NEWBUFF),NXMAX
                  WRITE(*,101) 'FATAL ERROR: NAXIS(1).GT.NXMAX'
                  CALL Deallocate_Array_IMAGEN
                  CALL Deallocate_Array_IMAGEN_
                  STOP
                END IF
                IF(NAXIS(2,NEWBUFF).GT.NYMAX)THEN
                  WRITE(*,100) 'NAXIS(2), NYMAX: '
                  WRITE(*,*) NAXIS(2,NEWBUFF),NYMAX
                  WRITE(*,101) 'FATAL ERROR: NAXIS(2).GT.NYMAX'
                  CALL Deallocate_Array_IMAGEN
                  CALL Deallocate_Array_IMAGEN_
                  STOP
                END IF
                DO I=1,NAXIS(2,NEWBUFF)
                  DO J=1,NAXIS(1,NEWBUFF)
                    IMAGEN(J,I,NEWBUFF)=0.
                  END DO
                END DO
                IF(CINERR.NE.'x')THEN
                  NAXIS(1,NEWBUFF+NMAXBUFF/2)=3*NSIZEFB9(NEWBUFF)
                  NAXIS(2,NEWBUFF+NMAXBUFF/2)=3*NSIZEFB9(NEWBUFF)
                  NFRAMES(NEWBUFF+NMAXBUFF/2)=NFRAMES(NEWBUFF)
                  DO I=1,NAXIS(2,NEWBUFF+NMAXBUFF/2)
                    DO J=1,NAXIS(1,NEWBUFF+NMAXBUFF/2)
                      IMAGEN(J,I,NEWBUFF+NMAXBUFF/2)=0.
                    END DO
                  END DO
                END IF
                OPEN(11,FILE=INFILE_,STATUS='OLD',FORM='FORMATTED')
                READ(11,*) NFRAMES(NEWBUFF)
                DO IBOX9=1,NFRAMES(NEWBUFF)
                  CALL BUTTON(32,'combine',0)
                  CALL BUTTON(33,'shiftb9',0)
                  READ(11,101) INFILEBOX9
                  WRITE(*,100) 'Reading file '
                  WRITE(*,100) INFILEBOX9(1:TRUELEN(INFILEBOX9))
                  WRITE(*,100) '...'
                  CALL SLEEFITS(INFILEBOX9,.FALSE.,0,0,LBOX9,0)
! Note: the pattern of the frames in box-9 is the following:
!       6 9 4
!       3 1 7
!       8 5 2
                  IF(IBOX9.EQ.1)THEN
!!!                 I1=256
!!!                 J1=256
                    I1=1
                    J1=1
                  ELSEIF(IBOX9.EQ.2)THEN
!!!                 I1=000
!!!                 J1=512
                    I1=0
                    J1=2
                  ELSEIF(IBOX9.EQ.3)THEN
!!!                 I1=256
!!!                 J1=000
                    I1=1
                    J1=0
                  ELSEIF(IBOX9.EQ.4)THEN
!!!                 I1=512
!!!                 J1=512
                    I1=2
                    J1=2
                  ELSEIF(IBOX9.EQ.5)THEN
!!!                 I1=000
!!!                 J1=256
                    I1=0
                    J1=1
                  ELSEIF(IBOX9.EQ.6)THEN
!!!                 I1=512
!!!                 J1=000
                    I1=2
                    J1=0
                  ELSEIF(IBOX9.EQ.7)THEN
!!!                 I1=256
!!!                 J1=512
                    I1=1
                    J1=2
                  ELSEIF(IBOX9.EQ.8)THEN
!!!                 I1=000
!!!                 J1=000
                    I1=0
                    J1=0
                  ELSEIF(IBOX9.EQ.9)THEN
!!!                 I1=512
!!!                 J1=256
                    I1=2
                    J1=1
                  END IF
                  I1=I1*NSIZEFB9(NEWBUFF)
                  J1=J1*NSIZEFB9(NEWBUFF)
                  DO I=1,NAXISFRAME(2,IBOX9,NEWBUFF)
                    DO J=1,NAXISFRAME(1,IBOX9,NEWBUFF)
                      IMAGEN(J+J1,I+I1,NEWBUFF)=IMAGEN_(J,I)
                    END DO
                  END DO
                  IF(CINERR.EQ.'g')THEN
                    DO I=1,NAXISFRAME(2,IBOX9,NEWBUFF)
                      DO J=1,NAXISFRAME(1,IBOX9,NEWBUFF)
                        IMAGEN(J+J1,I+I1,NEWBUFF+NMAXBUFF/2)= &
                         SQRT(ABS(IMAGEN(J+J1,I+I1,NEWBUFF))/GAIN_ARRAY+REAL(NCOADDS_ARRAY)*RNOISE_ARRAY*RNOISE_ARRAY)
                      END DO
                    END DO
                  ELSEIF(CINERR.EQ.'f')THEN
                    CALL GUESSEF(INFILEBOX9,ERRFILE)
                    WRITE(*,100) ERRFILE(1:TRUELEN(ERRFILE))
                    WRITE(*,100) '...'
                    CALL SLEEFITS(ERRFILE,.FALSE.,0,0,LBOX9,0)
                    DO I=1,NAXISFRAME(2,IBOX9,NEWBUFF)
                      DO J=1,NAXISFRAME(1,IBOX9,NEWBUFF)
                        IMAGEN(J+J1,I+I1,NEWBUFF+NMAXBUFF/2)=IMAGEN_(J,I)
                      END DO
                    END DO
                  END IF
                  WRITE(*,101) '   ...OK!'
                END DO
                CLOSE(11)
              ELSE
                WRITE(*,100) 'Reading file: '
                WRITE(*,101) INFILE_
                CALL SLEEFITS(INFILE_,.FALSE.,IROTATE,NEWBUFF,LBOX9,NEWBUFF)
                WRITE(*,101) ' ...OK!'
                IF(ANYNULL) THEN
                  WRITE(*,101) '***WARNING***'
                  WRITE(*,101) '=> ANYNULL=.TRUE.!'
                  !CALL PGEND
                  !STOP
                END IF
                NFRAMES(NEWBUFF)=1
                NSIZEFB9(NEWBUFF)=0
                NAXISFRAME(1,1,NEWBUFF)=0
                NAXISFRAME(2,1,NEWBUFF)=0
              END IF
              NCBUFF=NEWBUFF
              INFILE(NCBUFF)=INFILE_
              DO K=1,NMAXBUFF
                IF(K.LE.NMAXBUFF/2)THEN
                  WRITE(CDUMMY,'(A2,I1,A1)') '#[',K,']'
                  CALL RMBLANK(CDUMMY,CDUMMY,L)
                  IF(K.EQ.NCBUFF)THEN
                    CALL BUTTON(34+K,CDUMMY(1:L),1)
                  ELSE
                    CALL BUTTON(34+K,CDUMMY(1:L),0)
                  END IF
                ELSE
                  WRITE(CDUMMY,'(A5,I1)') 'err #',K-NMAXBUFF/2
                  L=TRUELEN(CDUMMY)
                  IF(K.EQ.NCBUFF)THEN
                    CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),1)
                  ELSE
                    CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),0)
                  END IF
                END IF
              END DO
              NX1=1
              NX2=NAXIS(1,NCBUFF)
              NY1=1
              NY2=NAXIS(2,NCBUFF)
              IF(LBOX9)THEN
                CALL STATISTICB9(NFRAMES(NCBUFF),NCBUFF,FMEAN,FSIGMA)
              ELSE
                CALL STATISTIC(NCBUFF,NX1,NX2,NY1,NY2,.FALSE.,.TRUE.,.TRUE.,0.0,.FALSE.)
              END IF
              IF(LFIRSTPLOT)THEN
                CALL ACTIVEBUT(JUST,MODECUT)
                IF(FSIGMA.GT.0.0)THEN
                  BG=FMEAN-5.*FSIGMA
                  FG=FMEAN+5.*FSIGMA
                ELSE
                  BG=FMEAN-1.0
                  FG=FMEAN+1.0
                END IF
                CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              ELSE
                CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
              END IF
              IF(LSTACK.AND.(NSTACK.GT.0))THEN
                CALL PGSCI(3)
                CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
                CALL PGSCI(1)
              END IF
              IF(LPBOUND) CALL PBOUND
              IF(LPMAPP) CALL PMAPP
              IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
              IF(LFIRSTPLOT)THEN
                CALL HISTOGRAM(NCBUFF)
                LFIRSTPLOT=.FALSE.
              END IF
              NY_CUT=NY2/2
              IF(NY_CUT.EQ.0) NY_CUT=1
              DO J=1,NAXIS(1,NCBUFF)
                XCUTX(J)=REAL(J)
                YCUTX(J)=IMAGEN(J,NY_CUT,NCBUFF)
              END DO
              WRITE(CDUMMY,'(A1,I10,A1,I10,A1)') '[',NY_CUT,',',NY_CUT,']'
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL SUBPLOT(NAXIS(1,NCBUFF),1,NAXIS(1,NCBUFF),XCUTX,YCUTX,XCUTX,YCUTX,.TRUE.,.TRUE.,.FALSE.,.FALSE., &
               'x axis','signal',CDUMMY(1:L),NCBUFF,201,1.0)
              IF(LOVERCUTS)THEN
                CALL BUTTON(51,'over=TRUE',0)
                CALL BUTTON(51,'over=TRUE',NCOLOR4)
              ELSE
                CALL BUTTON(51,'over=FALSE',0)
                CALL BUTTON(51,'over=FALSE',NCOLOR4)
              END IF
!..............................................................................
              IF((CINERR.EQ.'x').OR.(LBOX9))THEN               !no hacemos nada
!..............................................................................
              ELSEIF(CINERR.EQ.'0')THEN
                NAXIS(1,NCBUFF+NMAXBUFF/2)=NAXIS(1,NCBUFF)
                NAXIS(2,NCBUFF+NMAXBUFF/2)=NAXIS(2,NCBUFF)
                NFRAMES(NCBUFF+NMAXBUFF/2)=NFRAMES(NCBUFF)
                DO I=1,NAXIS(2,NCBUFF+NMAXBUFF/2)
                  DO J=1,NAXIS(1,NCBUFF+NMAXBUFF/2)
                    IMAGEN(J,I,NCBUFF+NMAXBUFF/2)=0.
                  END DO
                END DO
!..............................................................................
! Para calcular los errores, es importante notar que el numero de coadds es un
! factor importante; de hecho, el ruido debido a la sen~al hay que calcularlo
! a partir de la sen~al en cada coadd, no a partir de la sen~al total integrada
! que, logicamente, es mucho mayor. Notar que en la formula que sigue, el
! factor NCOADSS_ARRAY que multiplica al error en cada coadd se cancela con
! el mismo factor que aparece al dividir la sen~al total para obtener la
! sen~al en cada coadd. Sin embargo, este factor no se cancela en el termino
! correspondiente al ruido de lectura.
              ELSEIF(CINERR.EQ.'g')THEN
                NAXIS(1,NCBUFF+NMAXBUFF/2)=NAXIS(1,NCBUFF)
                NAXIS(2,NCBUFF+NMAXBUFF/2)=NAXIS(2,NCBUFF)
                NFRAMES(NCBUFF+NMAXBUFF/2)=NFRAMES(NCBUFF)
                DO I=1,NAXIS(2,NCBUFF+NMAXBUFF/2)
                  DO J=1,NAXIS(1,NCBUFF+NMAXBUFF/2)
                    IMAGEN(J,I,NCBUFF+NMAXBUFF/2)= &
                     SQRT(ABS(IMAGEN(J,I,NCBUFF))/GAIN_ARRAY+REAL(NCOADDS_ARRAY)*RNOISE_ARRAY*RNOISE_ARRAY)
                  END DO
                END DO
!..............................................................................
              ELSEIF(CINERR.EQ.'f')THEN
                LOGFILE=.FALSE.
                CALL GUESSEF(INFILE_,ERRFILE)
                DO WHILE(.NOT.LOGFILE)
                  ERRFILE=READC('Input FITS error file name',ERRFILE,'@')
                  IF((INDEX(ERRFILE,'*').NE.0).OR.(INDEX(ERRFILE,'?').NE.0))THEN
                    L1=TRUEBEG(ERRFILE)
                    L2=TRUELEN(ERRFILE)
                    ISYSTEM=SYSTEMFUNCTION('ls '//ERRFILE(L1:L2))
                  ELSE
                    INQUIRE(FILE=ERRFILE,EXIST=LOGFILE)
                    IF(.NOT.LOGFILE)THEN
                      WRITE(*,101) '***ERROR***'
                      WRITE(*,100) '=> This file does not exist.'
                      WRITE(*,100) ' Try again. (press <CR>...)'
                      READ(*,*)
                    END IF
                  END IF
                END DO
                LASK=.TRUE.
                DO WHILE(LASK)
                  WRITE(CDUMMY,*) IROTATE
                  IROTATE=READI('Rotation angle to be applied',CDUMMY)
                  LASK=((IROTATE.NE.0).AND.(IABS(IROTATE).NE.90))
                  IF(LASK)THEN
                    WRITE(*,101) '***ERROR***'
                    WRITE(*,100) '=> Invalid angle.'
                    WRITE(*,101) ' Only 0, +/-90, +/-180 +/-270'
                    WRITE(*,101) 'deg. are valid options.'
                    WRITE(*,100) '(press <CR> to continue...)'
                    READ(*,*)
                  END IF
                END DO
                WRITE(*,100) 'Reading file: '
                WRITE(*,101) ERRFILE
                CALL SLEEFITS(ERRFILE,.FALSE.,IROTATE,NCBUFF+NMAXBUFF/2,LBOX9,0)
                WRITE(*,101) ' ...OK!'
                IF(ANYNULL) THEN
                  WRITE(*,101) '***WARNING***'
                  WRITE(*,100) '=> ANYNULL=.TRUE.!'
                  !CALL PGEND
                  !STOP
                END IF
                INFILE(NCBUFF+NMAXBUFF/2)=INFILE_
              END IF
!..............................................................................
            END IF
!
            CALL BUTTON(NB,'[l]oad',0)
            CALL BUTTON(NB,'[l]oad',NCOLOR1)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.2)THEN
            CALL BUTTON(NB,'[s]ave',5)
            CSAVE(1:1)=READC('[f/F]its or [r]educeme format (f/F/r/x)',CSAVE,'fFrx')
            IF(CSAVE.NE.'x')THEN
              IF(NCBUFF.LE.NMAXBUFF/2)THEN
                CSAVEERR(1:1)=READC('Save error frame (y/n)',CSAVEERR,'yn')
              ELSE
                CSAVEERR='n'
              END IF
              LOGFILE=.TRUE.
              LOGFILERR=.FALSE.
              DO WHILE(LOGFILE)
                IF((CSAVE.EQ.'f').OR.(CSAVE.EQ.'F'))THEN
                  OUTFILE=READC('Output FITS file name (none=EXIT)','*.fits','@')
                ELSE
                  OUTFILE=READC('Output REDUCEME file name (none=EXIT)','*.u','@')
                END IF
                IF((INDEX(OUTFILE,'*').NE.0).OR.(INDEX(OUTFILE,'?').NE.0))THEN
                  L1=TRUEBEG(OUTFILE)
                  L2=TRUELEN(OUTFILE)
                  ISYSTEM=SYSTEMFUNCTION('ls '//OUTFILE(L1:L2))
                ELSEIF(OUTFILE.EQ.'none')THEN
                  LOGFILE=.FALSE.
                  LOGFILERR=.TRUE.
                ELSE
                  INQUIRE(FILE=OUTFILE,EXIST=LOGFILE)
                  IF(LOGFILE)THEN
                    WRITE(*,101) '***ERROR***'
                    WRITE(*,100) '=> This file does already exist.'
                    WRITE(*,100) ' Try again. (press <CR>...)'
                    READ(*,*)
                  END IF
                END IF
              END DO
              IF(.NOT.LOGFILERR)THEN
                IF(CSAVE.EQ.'F')THEN
                  LOGFILE=.FALSE.
                  LOGFILERR=.FALSE.
                  DO WHILE(.NOT.LOGFILE)
                    WRITE(*,101)'* Including header from other file'
                    HDRFILE=READC('Header FITS file name (none=EXIT)','*.fits','@')
                    IF((INDEX(HDRFILE,'*').NE.0).OR.(INDEX(HDRFILE,'?').NE.0))THEN
                      L1=TRUEBEG(HDRFILE)
                      L2=TRUELEN(HDRFILE)
                      ISYSTEM=SYSTEMFUNCTION('ls '//HDRFILE(L1:L2))
                    ELSEIF(HDRFILE.EQ.'none')THEN
                      LOGFILE=.TRUE.
                      LOGFILERR=.TRUE.
                    ELSE
                      INQUIRE(FILE=HDRFILE,EXIST=LOGFILE)
                      IF(.NOT.LOGFILE)THEN
                        WRITE(*,101) '***ERROR***'
                        WRITE(*,100) '=> This file does not exist.'
                        WRITE(*,100) ' Try again. (press <CR>...)'
                        READ(*,*)
                      END IF
                    END IF
                  END DO
                END IF
              END IF
              IF(.NOT.LOGFILERR)THEN               !the file name is not 'none'
                WRITE(*,101) 'Writting file...'
                IF(CSAVE.EQ.'f')THEN
                  CALL SESCRFITS(OUTFILE,NCBUFF,'none')
                ELSEIF(CSAVE.EQ.'F')THEN
                  CALL SESCRFITS(OUTFILE,NCBUFF,HDRFILE)
                ELSE
                  CALL SESCRREDU(OUTFILE,NCBUFF,.FALSE.)
                END IF
                WRITE(*,101) ' ...OK!'
                IF(CSAVEERR.EQ.'y')THEN
                  LOGFILE=.TRUE.
                  CALL GUESSEF(OUTFILE,ERRFILE)
                  DO WHILE(LOGFILE)
                    IF((CSAVE.EQ.'f').OR.(CSAVE.EQ.'F'))THEN
                      ERRFILE=READC('Output FITS error file name',ERRFILE,'@')
                    ELSE
                      ERRFILE=READC('Output REDUCEME error file name',ERRFILE,'@')
                    END IF
                    IF((INDEX(ERRFILE,'*').NE.0).OR.(INDEX(ERRFILE,'?').NE.0))THEN
                      L1=TRUEBEG(ERRFILE)
                      L2=TRUELEN(ERRFILE)
                      ISYSTEM=SYSTEMFUNCTION('ls '//ERRFILE(L1:L2))
                    ELSE
                      INQUIRE(FILE=ERRFILE,EXIST=LOGFILE)
                      IF(LOGFILE)THEN
                        WRITE(*,101) '***ERROR***'
                        WRITE(*,100) '=> This file does already exist.'
                        WRITE(*,100) ' Try again. (press <CR>...)'
                        READ(*,*)
                      END IF
                    END IF
                  END DO
                  WRITE(*,100) 'Writting file...'
                  IF(CSAVE.EQ.'f')THEN
                    CALL SESCRFITS(ERRFILE,NCBUFF+NMAXBUFF/2,'none')
                  ELSEIF(CSAVE.EQ.'F')THEN
                    CALL SESCRFITS(ERRFILE,NCBUFF+NMAXBUFF/2,HDRFILE)
                  ELSE
                    CALL SESCRREDU(ERRFILE,NCBUFF+NMAXBUFF/2,.TRUE.)
                  END IF
                  WRITE(*,101) ' ...OK!'
                END IF
              END IF
            END IF
            CALL BUTTON(NB,'[s]ave',0)
            CALL BUTTON(NB,'[s]ave',NCOLOR1)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.3)THEN
            CALL BUTTON(NB,'[b]oundary',5)
            CALL BOUNDARY(LBOUNDARY)
            IF(LBOUNDARY)THEN
              CALL BUTTON(4,'find lines 1',NCOLOR5)
              CALL BUTTON(14,'find lines 2',NCOLOR5)
              LPBOUND=.TRUE.
              CALL BUTTON(13,'plot bou.',0)
              CALL BUTTON(13,'plot bou.',1)
            END IF
            CALL BUTTON(NB,'[b]oundary',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.4)THEN
            CALL BUTTON(NB,'find lines 1',5)
            WRITE(*,101) 'Press mouse within boundary...'
            CALL PGSCI(5)
            CALL RPGBAND(7,0,0.,0.,XC,YC,CH)
            CALL PGSCI(1)
            IF(CH.EQ.'X')THEN
              WRITE(*,101) '...operation canceled!'
            ELSE
              WRITE(*,101) '...OK!'
              CALL WHEREAMI(XC,YC,COEFFBL00,COEFFBA00,-1,-1,.TRUE.)
              WRITE(CDUMMY,*) NYFINDAL
              NYFINDAL=READILIM('Y grid looking for spectral lines',CDUMMY,1,NXYMAX/4)
              WRITE(CDUMMY,*) NINLINE
              NINLINE=READILIM('No. of peaks/line found',CDUMMY,1,NYFINDAL)
              WRITE(CDUMMY,*) NWIDTH
              NWIDTH=READILIM('No. of pixels to search for peaks (odd)',CDUMMY,3,NWIDTHMAX)
!!!           CALL EXPLINE(XC,YC,COEFFBL00,COEFFBA00)
              CALL EXPLINE(XC,COEFFBL00,COEFFBA00)
            END IF
            CALL BUTTON(NB,'find lines 1',0)
            CALL BUTTON(NB,'find lines 1',NCOLOR5)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.5)THEN
            CALL BUTTON(NB,'[z]oom',5)
            CALL SUBZOOM(NX1,NX2,NY1,NY2,LOK)
            IF(LOK)THEN
              CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              IF(LSTACK.AND.(NSTACK.GT.0))THEN
                CALL PGSCI(3)
                CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
                CALL PGSCI(1)
              END IF
              IF(LPBOUND) CALL PBOUND
              IF(LPMAPP) CALL PMAPP
              IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
!!!              CALL STATISTIC(NCBUFF,NX1,NX2,NY1,NY2,
!!!               .FALSE.,.TRUE.,.TRUE.,0.0,.FALSE.)
            END IF
            CALL BUTTON(NB,'[z]oom',0)
            CALL BUTTON(NB,'[z]oom',NCOLOR2)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.6)THEN
            CALL BUTTON(NB,'[w]hole',5)
            IF((NAXIS(1,NCBUFF).EQ.0).OR.(NAXIS(2,NCBUFF).EQ.0))THEN
              WRITE(*,101) '***ERROR***'
              WRITE(*,100) '=> NCBUFF: '
              WRITE(*,*) NCBUFF
              WRITE(*,101) '=> This buffer has not been defined!'
            ELSE
              NX1=1
              NX2=NAXIS(1,NCBUFF)
              NY1=1
              NY2=NAXIS(2,NCBUFF)
              CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              IF(LSTACK.AND.(NSTACK.GT.0))THEN
                CALL PGSCI(3)
                CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
                CALL PGSCI(1)
              END IF
              IF(LPBOUND) CALL PBOUND
              IF(LPMAPP) CALL PMAPP
              IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            END IF
!!!            CALL STATISTIC(NCBUFF,NX1,NX2,NY1,NY2,
!!!             .FALSE.,.TRUE.,.TRUE.,0.0,.FALSE.)
            CALL BUTTON(NB,'[w]hole',0)
            CALL BUTTON(NB,'[w]hole',NCOLOR2)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.7)THEN
            IF(JUST.EQ.1)THEN
              JUST=0
              CALL BUTTON(NB,'  x![=]y',0)
              CALL BUTTON(NB,'  x![=]y',NCOLOR2)
            ELSE
              JUST=1
              CALL BUTTON(NB,'x[=]y',0)
              CALL BUTTON(NB,'x[=]y',NCOLOR2)
            END IF
            CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
            IF(LSTACK.AND.(NSTACK.GT.0))THEN
              CALL PGSCI(3)
              CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
              CALL PGSCI(1)
            END IF
            IF(LPBOUND) CALL PBOUND
            IF(LPMAPP) CALL PMAPP
            IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.8)THEN
            CALL BUTTON(NB,'[m]apping',5)
            CALL MAPPING(LBOUNDARY,LMAPPING)
            IF(LMAPPING)THEN
              IF(.NOT.LPMAPP) CALL BUTTON(18,'[p]lot mapp.',0)
              CALL BUTTON(9,'wa[v]el.cal.',NCOLOR5)
              CALL BUTTON(25,'[r]ectify',NCOLOR6)
              CALL BUTTON(26,'s[k]y sub.',NCOLOR6)
              CALL BUTTON(27,'X-[e]xtract',NCOLOR6)
              CALL BUTTON(28,'[u]nrectify',NCOLOR6)
            END IF
            CALL BUTTON(NB,'[m]apping',0)
            CALL BUTTON(NB,'[m]apping',NCOLOR5)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.9)THEN
            CALL BUTTON(NB,'wa[v]el.cal.',5)
            LOGFILE=.FALSE.
            LOGFILERR=.FALSE.
            DO WHILE(.NOT.LOGFILE)
              INFILE_=READC('Input file name with polynomial wavelength calibration (none=EXIT)','polw*','@')
              IF((INDEX(INFILE_,'*').NE.0).OR.(INDEX(INFILE_,'?').NE.0))THEN
                L1=TRUEBEG(INFILE_)
                L2=TRUELEN(INFILE_)
                ISYSTEM=SYSTEMFUNCTION('ls '//INFILE_(L1:L2))
              ELSEIF(INFILE_.EQ.'none')THEN
                LOGFILE=.TRUE.
                LOGFILERR=.TRUE.
              ELSE
                INQUIRE(FILE=INFILE_,EXIST=LOGFILE)
                IF(.NOT.LOGFILE)THEN
                  WRITE(*,101) '***ERROR***'
                  WRITE(*,100) '=> This file does not exist.'
                  WRITE(*,100) ' Try again. (press <CR>...)'
                  READ(*,*)
                END IF
              END IF
            END DO
!
            IF(.NOT.LOGFILERR)THEN                 !the file name is not 'none'
              CALL READPWAVE(INFILE_,NDEGWV,CWV)
              CALL BUTTON(19,'plot wave.',1)
              LPWAVE=.TRUE.
              WAVEFILE=INFILE_
            END IF
            CALL BUTTON(NB,'wa[v]el.cal.',0)
            CALL BUTTON(NB,'wa[v]el.cal.',NCOLOR5)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.10)THEN
            IF(LSTACK)THEN
              IF(NSTACK.GT.0)THEN
                WRITE(CDUMMY,*) NSTACK
                CALL RMBLANK(CDUMMY,CDUMMY,L1)
                CALL BUTTON(10,'stack '//CDUMMY(1:L1),5)
                WRITE(*,101) '(p) pause stack'
                WRITE(*,101) '(d) delete stack'
                WRITE(*,101) '(r) delete & restart stack'
                WRITE(*,101) '(1) fit y(x)'
                WRITE(*,101) '(2) fit x(y)'
                WRITE(*,101) '(3) define BL boundary'
                WRITE(*,101) '(4) define BA boundary'
                WRITE(*,101) '(0) exit'
                CDUM(1:1)=READC('Option','p','pdr01234')
!
                IF(CDUM.EQ.'p')THEN
                  LSTACK=.FALSE.
                  WRITE(CDUMMY,*) NSTACK
                  CALL RMBLANK(CDUMMY,CDUMMY,L1)
                  CALL BUTTON(10,'stack '//CDUMMY(1:L1),0)
!..............................................................................
                ELSEIF(CDUM.EQ.'d')THEN
                  LSTACK=.FALSE.
                  NSTACK=0
!..............................................................................
                  CALL BUTTON(10,'stack 0',0)
                ELSEIF(CDUM.EQ.'r')THEN
                  NSTACK=0
                  CALL BUTTON(10,'stack 0',0)
                  CALL BUTTON(10,'stack 0',1)
!..............................................................................
                ELSEIF((CDUM.EQ.'1').OR.(CDUM.EQ.'2'))THEN
                  NDEG=READILIM('Polynomial degree','0',0,MIN0(19,NSTACK-1))
                  IF(CDUM.EQ.'1')THEN
                    CALL DRAWPOLYX(NSTACK,XSTACK,YSTACK,NDEG,COEFF)
                  ELSE
                    CALL DRAWPOLXY(NSTACK,XSTACK,YSTACK,NDEG,COEFF)
                  END IF
                  WRITE(CDUMMY,*) NSTACK
                  CALL RMBLANK(CDUMMY,CDUMMY,L1)
                  CALL BUTTON(10,'stack '//CDUMMY(1:L1),0)
                  CALL BUTTON(10,'stack '//CDUMMY(1:L1),1)
!..............................................................................
                ELSEIF(CDUM.EQ.'3')THEN
                  NDEG=READILIM('Polynomial degree','0',0,MIN0(19,NSTACK-1))
                  CALL DRAWPOLYX(NSTACK,XSTACK,YSTACK,NDEG,COEFF)
                  CREFINE(1:1)=READC('Refine line fit (y/n)','y','yn')
                  IF(CREFINE.EQ.'y')THEN
                    CALL FINDMML(NSTACK,1,NSTACK,XSTACK,XMINFIT,XMAXFIT)
                    WRITE(CDUMMY,*) NWIDTH
                    NWIDTH=READILIM('No. of pixels to search for peaks (odd)',CDUMMY,3,NWIDTHMAX)
                    CALL REFINEBL(NDEG,COEFF,XMINFIT,XMAXFIT)
                    WRITE(*,101) '=> fit result: '
                    DO K=1,NDEG+1
                      WRITE(*,'(A2,I2,A3,$)') 'a(',K-1,'): '
                      WRITE(*,*) COEFF(K)
                    END DO
                  END IF
                  CUPDATE(1:1)=READC('Update BL lines (y/n)','y','yn')
                  IF(CUPDATE.EQ.'y')THEN
                    NLINBL=NLINBL+1
                    NDEGBL(NLINBL)=NDEG
                    DO K=1,NDEG+1
                      COEFFBL(K,NLINBL)=COEFF(K)
                    END DO
                    FSHIFT=READF('Shift in Y (+ upwards, - downwards)','0.0')
                    COEFFBL(1,NLINBL)=COEFFBL(1,NLINBL)+FSHIFT
                    LBOUNDL=(NLINBL.GE.2)
                    XMINBL(NLINBL)=0.
                    YMINBL(NLINBL)=FPOLY(NDEGBL(NLINBL),COEFFBL(1,NLINBL),XMINBL(NLINBL))
                    XMAXBL(NLINBL)=REAL(NAXIS(1,NCBUFF)+1)
                    YMAXBL(NLINBL)=FPOLY(NDEGBL(NLINBL),COEFFBL(1,NLINBL),XMAXBL(NLINBL))
                    CALL SORTLINES
                    WRITE(*,101) '***WARNING***'
                    WRITE(*,101) '=> RESET shrinking/stretching'
                    WRITE(*,*)
                  END IF
!..............................................................................
                ELSEIF(CDUM.EQ.'4')THEN
                  NDEG=READILIM('Polynomial degree','0',0,MIN0(19,NSTACK-1))
                  CALL DRAWPOLXY(NSTACK,XSTACK,YSTACK,NDEG,COEFF)
                  CREFINE(1:1)=READC('Refine line fit (y/n)','y','yn')
                  IF(CREFINE.EQ.'y')THEN
                    WRITE(CDUMMY,*) NWIDTH
                    NWIDTH=READILIM('No. of pixels to search for peaks (odd)',CDUMMY,3,NWIDTHMAX)
                    CALL REFINEBA(NDEG,COEFF,XSTACK(NSTACK/2),FMEAN_,FSIGMA_)
                    WRITE(*,101) '=> fit result: '
                    DO K=1,NDEG+1
                      WRITE(*,'(A2,I2,A3,$)') 'a(',K-1,'): '
                      WRITE(*,*) COEFF(K)
                    END DO
                  END IF
                  CUPDATE(1:1)=READC('Update BA lines (y/n)','y','yn')
                  IF(CUPDATE.EQ.'y')THEN
                    NLINBA=NLINBA+1
                    NDEGBA(NLINBA)=NDEG
                    DO K=1,NDEG+1
                      COEFFBA(K,NLINBA)=COEFF(K)
                    END DO
                    LBOUNDA=(NLINBA.GE.2)
                    CALL SORTLINES
                    WRITE(*,101) '***WARNING***'
                    WRITE(*,101) '=> RESET shrinking/stretching'
                    WRITE(*,*)
                  END IF
!..............................................................................
                ELSE
                  CALL BUTTON(10,'stack '//CDUMMY(1:L1),0)
                  CALL BUTTON(10,'stack '//CDUMMY(1:L1),1)
                END IF
!
                IF((CDUM.EQ.'3').OR.(CDUM.EQ.'4'))THEN
                  IF(LBOUNDL.OR.LBOUNDA)THEN
                    LPBOUND=.TRUE.
                    CALL BUTTON(13,'plot bou.',0)
                    CALL BUTTON(13,'plot bou.',1)
                  END IF
                  LBOUNDARY=(LBOUNDL.AND.LBOUNDA)
                  IF(LBOUNDARY)THEN
                    CALL BUTTON(4,'find lines 1',NCOLOR5)
                    CALL BUTTON(14,'find lines 2',NCOLOR5)
                  END IF
                  WRITE(CDUMMY,*) NSTACK
                  CALL RMBLANK(CDUMMY,CDUMMY,L1)
                  CALL BUTTON(10,'stack '//CDUMMY(1:L1),0)
                  CALL BUTTON(10,'stack '//CDUMMY(1:L1),1)
                END IF
              ELSE
                LSTACK=.FALSE.
                WRITE(CDUMMY,*) NSTACK
                CALL RMBLANK(CDUMMY,CDUMMY,L1)
                CALL BUTTON(10,'stack '//CDUMMY(1:L1),0)
              END IF
            ELSE
              LSTACK=.TRUE.
              WRITE(CDUMMY,*) NSTACK
              CALL RMBLANK(CDUMMY,CDUMMY,L1)
              CALL BUTTON(10,'stack '//CDUMMY(1:L1),0)
              CALL BUTTON(10,'stack '//CDUMMY(1:L1),1)
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.11)THEN
            CALL BUTTON(NB,'[q]uit',5)
            IF(CH.EQ.'Q')THEN
              LEXIT=.TRUE.
            ELSE
              CSURE(1:1)=READC('Do you really want to exit (y/n)','@','yn')
              IF(CSURE.EQ.'y')THEN
                LEXIT=.TRUE.
              ELSE
                CALL BUTTON(NB,'[q]uit',0)
                CALL BUTTON(NB,'[q]uit',NCOLOR1)
              END IF
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.12)THEN
            CALL BUTTON(NB,'[i]math',5)
            CALL SUBIMATH
            CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
            CALL BUTTON(NB,'[i]math',0)
            CALL BUTTON(NB,'[i]math',NCOLOR1)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.13)THEN
            IF(LPBOUND)THEN
              CALL BUTTON(13,'plot bou.',0)
              LPBOUND=.FALSE.
            ELSE
              CALL BUTTON(13,'plot bou.',1)
              LPBOUND=.TRUE.
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.14)THEN
            CALL BUTTON(NB,'find lines 2',5)
            WRITE(CDUMMY,*) NWIDTH
            NWIDTH=READILIM('No. of pixels to search for peaks (odd)',CDUMMY,3,NWIDTHMAX)
            CH='A'
            DO WHILE(CH.NE.'X')
              WRITE(*,101) 'Press mouse in sky line...'
              CALL PGSCI(5)
              CALL RPGBAND(7,0,0.,0.,XC,YC,CH)
              CALL PGSCI(1)
              IF(CH.EQ.'X')THEN
                WRITE(*,101) '...operation canceled!'
              ELSE
                CALL WHEREAMI(XC,YC,COEFFBL00,COEFFBA00,-1,-1,.FALSE.)
                CALL INTERSEC(NDEGBA00,COEFFBA00,NDEGBL(1),COEFFBL(1,1),XC,XMIN,YMIN)
                CALL INTERSEC(NDEGBA00,COEFFBA00,NDEGBL(NLINBL),COEFFBL(1,NLINBL),XC,XMAX,YMAX)
                DO I=1,NXYMAX
                  YP(I)=YMIN+REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)
                  XP(I)=FPOLY(NDEGBA00,COEFFBA00,YP(I))
                END DO
                CALL PGSCI(7)
                CALL PGLINE(NXYMAX,XP,YP)
                CALL PGSCI(1)
                CALL REFINEBA(NDEGBA00,COEFFBA00,XP(NXYMAX/2),FMEAN_,FSIGMA_)
                WRITE(*,100) 'Press left mouse button to add new line'
                WRITE(*,101) ' to the pool of BA lines'
                CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
                IF(CH.EQ.'A')THEN
                  NLINBA=NLINBA+1
                  NDEGBA(NLINBA)=NDEGBA00
                  DO K=1,NDEGBA00+1
                    COEFFBA(K,NLINBA)=COEFFBA00(K)
                  END DO
                  CALL SORTLINES
                  WRITE(*,101) 'Last line has been added to the pool!'
                  WRITE(*,100) '* New total number of BA lines: '
                  WRITE(*,*) NLINBA
                  WRITE(*,*)
                ELSE
                  CH='A'
                END IF
              END IF
            END DO
            CALL BUTTON(NB,'find lines 2',0)
            CALL BUTTON(NB,'find lines 2',NCOLOR5)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.15)THEN
            IF(MODECUT.EQ.1)THEN
              CALL BUTTON(NB,'single',-1)
              CALL BUTTON(NB,'region',0)
              CALL BUTTON(NB,'region',NCOLOR2)
              MODECUT=2
            ELSEIF(MODECUT.EQ.2)THEN
              CALL BUTTON(NB,'region',-1)
              CALL BUTTON(NB,'all',0)
              CALL BUTTON(NB,'all',NCOLOR2)
              MODECUT=3
            ELSE
              CALL BUTTON(NB,'all',-1)
              CALL BUTTON(NB,'single',0)
              CALL BUTTON(NB,'single',NCOLOR2)
              MODECUT=1
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.16)THEN
            CALL BUTTON(NB,'[x] cut',5)
            CALL PCUT('X',MODECUT,K1,K2,XP,YP)
            CALL BUTTON(NB,'[x] cut',0)
            CALL BUTTON(NB,'[x] cut',NCOLOR2)
            IF(LOVERCUTS)THEN
              CALL BUTTON(51,'over=TRUE',0)
              CALL BUTTON(51,'over=TRUE',NCOLOR4)
            ELSE
              CALL BUTTON(51,'over=FALSE',0)
              CALL BUTTON(51,'over=FALSE',NCOLOR4)
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.17)THEN
            CALL BUTTON(NB,'[y] cut',5)
            CALL PCUT('Y',MODECUT,K1,K2,XP,YP)
            CALL BUTTON(NB,'[y] cut',0)
            CALL BUTTON(NB,'[y] cut',NCOLOR2)
            IF(LOVERCUTS)THEN
              CALL BUTTON(51,'over=TRUE',0)
              CALL BUTTON(51,'over=TRUE',NCOLOR4)
            ELSE
              CALL BUTTON(51,'over=FALSE',0)
              CALL BUTTON(51,'over=FALSE',NCOLOR4)
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.18)THEN
            IF(LPMAPP)THEN
              CALL BUTTON(NB,'[p]lot mapp.',0)
              LPMAPP=.FALSE.
            ELSE
              CALL BUTTON(NB,'[p]lot mapp.',1)
              LPMAPP=.TRUE.
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.19)THEN
            IF(LPWAVE)THEN
              CALL BUTTON(NB,'plot wave.',0)
              LPWAVE=.FALSE.
            ELSE
              CALL BUTTON(NB,'plot wave.',1)
              LPWAVE=.TRUE.
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.20)THEN
            CALL BUTTON(NB,'s[Tt]atistic',5)
            IF(CH.NE.'T')THEN
              WRITE(CDUMMY,*) NBOX
              NBOX=READILIM('NBOX (-4=mask,-3=box9, -2=quad., -1=mouse, 0=keyboard)',CDUMMY,-4,MIN(NXMAX,NYMAX))
            END IF
            IF(NBOX.EQ.-4)THEN
              WRITE(CDUMMY,*) MASKVALUE
              MASKVALUE=READF('Image value to be ignored in statistics',CDUMMY)
              WRITE(CDUMMY,'(I10,A1,I10)') NX1,',',NX2
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL READ2I('Channel region ',CDUMMY(1:L),NNX1,NNX2)
              WRITE(CDUMMY,'(I10,A1,I10)') NY1,',',NY2
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL READ2I('Scan region ',CDUMMY(1:L),NNY1,NNY2)
              IF((NNX1.GE.1).AND.(NNX2.LE.NAXIS(1,NCBUFF)).AND.(NNY1.GE.1).AND.(NNY2.LE.NAXIS(2,NCBUFF)))THEN
                CALL PGSCI(3)
                CALL PGMOVE(REAL(NNX1),REAL(NNY1))
                CALL PGDRAW(REAL(NNX2),REAL(NNY1))
                CALL PGDRAW(REAL(NNX2),REAL(NNY2))
                CALL PGDRAW(REAL(NNX1),REAL(NNY2))
                CALL PGDRAW(REAL(NNX1),REAL(NNY1))
                CALL PGSCI(1)
                CALL STATISTIC(NCBUFF,NNX1,NNX2,NNY1,NNY2,.TRUE.,.TRUE.,.TRUE.,MASKVALUE,.TRUE.)
              END IF
            ELSEIF(NBOX.EQ.-3)THEN
              DO I=1,3
                IF(I.EQ.1)THEN
                  NNY1_=513
                ELSEIF(I.EQ.2)THEN
                  NNY1_=257
                ELSEIF(I.EQ.3)THEN
                  NNY1_=1
                END IF
                DO J=1,3
                  IF(J.EQ.1)THEN
                    NNX1_=1
                  ELSEIF(J.EQ.2)THEN
                    NNX1_=257
                  ELSEIF(J.EQ.3)THEN
                    NNX1_=513
                  END IF
                  DO K=1,4
                    IF(K.EQ.1)THEN
                      NNY1=NNY1_+128
                      NNX1=NNX1_
                    ELSEIF(K.EQ.2)THEN
                      NNY1=NNY1_+128
                      NNX1=NNX1_+128
                    ELSEIF(K.EQ.3)THEN
                      NNY1=NNY1_
                      NNX1=NNX1_
                    ELSEIF(K.EQ.4)THEN
                      NNY1=NNY1_
                      NNX1=NNX1_+128
                    END IF
                    NNY2=NNY1+127
                    NNX2=NNX1+127
                    CALL PGSCI(K+1)
                    CALL PGMOVE(REAL(NNX1),REAL(NNY1))
                    CALL PGDRAW(REAL(NNX2),REAL(NNY1))
                    CALL PGDRAW(REAL(NNX2),REAL(NNY2))
                    CALL PGDRAW(REAL(NNX1),REAL(NNY2))
                    CALL PGDRAW(REAL(NNX1),REAL(NNY1))
                    CALL PGSCI(1)
                    CALL STATISTIC(NCBUFF,NNX1,NNX2,NNY1,NNY2,.TRUE.,.TRUE.,.TRUE.,0.0,.FALSE.)
                    WRITE(*,100) 'Press <CR> to continue...'
                    READ(*,*)
                  END DO
                END DO
              END DO
            ELSEIF(NBOX.EQ.-2)THEN
              DO K=1,4
                IF(K.EQ.1)THEN
                  NNY1=NAXIS(2,NCBUFF)/2
                  NNX1=1
                ELSEIF(K.EQ.2)THEN
                  NNY1=NAXIS(2,NCBUFF)/2
                  NNX1=NAXIS(1,NCBUFF)/2
                ELSEIF(K.EQ.3)THEN
                  NNY1=1
                  NNX1=1
                ELSEIF(K.EQ.4)THEN
                  NNY1=1
                  NNX1=NAXIS(1,NCBUFF)/2
                END IF
                NNY2=NNY1+NAXIS(2,NCBUFF)/2-1
                NNX2=NNX1+NAXIS(1,NCBUFF)/2-1
                CALL PGSCI(K+1)
                CALL PGMOVE(REAL(NNX1),REAL(NNY1))
                CALL PGDRAW(REAL(NNX2),REAL(NNY1))
                CALL PGDRAW(REAL(NNX2),REAL(NNY2))
                CALL PGDRAW(REAL(NNX1),REAL(NNY2))
                CALL PGDRAW(REAL(NNX1),REAL(NNY1))
                CALL PGSCI(1)
                CALL STATISTIC(NCBUFF,NNX1,NNX2,NNY1,NNY2,.TRUE.,.TRUE.,.TRUE.,0.0,.FALSE.)
                WRITE(*,100) 'Press <CR> to continue...'
                READ(*,*)
              END DO
            ELSEIF(NBOX.EQ.-1)THEN
              NNX1=NX1
              NNX2=NX2
              NNY1=NY1
              NNY2=NY2
              CALL SUBZOOM(NNX1,NNX2,NNY1,NNY2,LOK)
              IF(LOK)THEN
                CALL PGSCI(3)
                CALL PGMOVE(REAL(NNX1),REAL(NNY1))
                CALL PGDRAW(REAL(NNX2),REAL(NNY1))
                CALL PGDRAW(REAL(NNX2),REAL(NNY2))
                CALL PGDRAW(REAL(NNX1),REAL(NNY2))
                CALL PGDRAW(REAL(NNX1),REAL(NNY1))
                CALL PGSCI(1)
                CALL STATISTIC(NCBUFF,NNX1,NNX2,NNY1,NNY2,.TRUE.,.TRUE.,.TRUE.,0.0,.FALSE.)
              END IF
            ELSEIF(NBOX.EQ.0)THEN
              WRITE(CDUMMY,'(I10,A1,I10)') NX1,',',NX2
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL READ2I('Channel region ',CDUMMY(1:L),NNX1,NNX2)
              WRITE(CDUMMY,'(I10,A1,I10)') NY1,',',NY2
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL READ2I('Scan region ',CDUMMY(1:L),NNY1,NNY2)
              IF((NNX1.GE.1).AND.(NNX2.LE.NAXIS(1,NCBUFF)).AND.(NNY1.GE.1).AND.(NNY2.LE.NAXIS(2,NCBUFF)))THEN
                CALL PGSCI(3)
                CALL PGMOVE(REAL(NNX1),REAL(NNY1))
                CALL PGDRAW(REAL(NNX2),REAL(NNY1))
                CALL PGDRAW(REAL(NNX2),REAL(NNY2))
                CALL PGDRAW(REAL(NNX1),REAL(NNY2))
                CALL PGDRAW(REAL(NNX1),REAL(NNY1))
                CALL PGSCI(1)
                CALL STATISTIC(NCBUFF,NNX1,NNX2,NNY1,NNY2,.TRUE.,.TRUE.,.TRUE.,0.0,.FALSE.)
              END IF
            ELSE
              CH='A'
              DO WHILE(CH.NE.'X')
                CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
                IF(CH.NE.'X')THEN
                  NNX1=NINT(XC)-NBOX/2
                  IF(NNX1.LT.1) NNX1=1
                  NNX2=NNX1+NBOX-1
                  IF(NNX2.GT.NAXIS(1,NCBUFF))THEN
                    NNX2=NAXIS(1,NCBUFF)
                    NNX1=NNX2-NBOX+1
                  END IF
                  NNY1=NINT(YC)-NBOX/2
                  IF(NNY1.LT.1) NNY1=1
                  NNY2=NNY1+NBOX-1
                  IF(NNY2.GT.NAXIS(2,NCBUFF))THEN
                    NNY2=NAXIS(2,NCBUFF)
                    NNY1=NNY2-NBOX+1
                  END IF
                  CALL PGSCI(3)
                  CALL PGMOVE(REAL(NNX1),REAL(NNY1))
                  CALL PGDRAW(REAL(NNX2),REAL(NNY1))
                  CALL PGDRAW(REAL(NNX2),REAL(NNY2))
                  CALL PGDRAW(REAL(NNX1),REAL(NNY2))
                  CALL PGDRAW(REAL(NNX1),REAL(NNY1))
                  CALL PGSCI(1)
                  CALL STATISTIC(NCBUFF,NNX1,NNX2,NNY1,NNY2,.TRUE.,.TRUE.,.TRUE.,0.0,.FALSE.)
                END IF
              END DO
            END IF
            CALL BUTTON(NB,'s[Tt]atistic',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.21)THEN
            CALL BUTTON(NB,'[?]buffers',5)
            DO K=1,NMAXBUFF
              IF((NAXIS(1,K).NE.0).AND.(NAXIS(2,K).NE.0))THEN
                IF(K.LE.NMAXBUFF/2)THEN
                  WRITE(*,'(A4,I1,A13,$)') '  #[',K,'] <DEFINED>: '
                ELSE
                  WRITE(*,'(A5,I1,A12,$)') 'err #',K-NMAXBUFF/2,' <DEFINED>: '
                END IF
                L=TRUELEN(INFILE(K))
                WRITE(*,100) INFILE(K)(1:L)
                WRITE(CDUMMY,'(A1,I6,A1,I6,A1)') '[',NAXIS(1,K),',',NAXIS(2,K),']'
                CALL RMBLANK(CDUMMY,CDUMMY,L)
                WRITE(*,101) ' --> '//CDUMMY(1:L)
              ELSE
                IF(K.LE.NMAXBUFF/2)THEN
                  WRITE(*,'(A4,I1,A13)') '  #[',K,'] <UNDEFINED>'
                ELSE
                  WRITE(*,'(A5,I1,A12)') 'err #',K-NMAXBUFF/2,' <UNDEFINED>'
                END IF
              END IF
            END DO
            CALL BUTTON(NB,'[?]buffers',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.22)THEN
            CALL BUTTON(NB,'postscript',5)
            NEWDEV=READC('New graphics device','/ps','@')
            IDNEW=PGOPEN(NEWDEV)
            IF(IDNEW.LE.0)THEN
              WRITE(*,101) 'ERROR: invalid graphics device.'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
              CALL PGSLCT(ID)
            ELSE
              WRITE(*,100) 'Creating new plot...'
              CALL SUBLOOK(.FALSE.,NCBUFF,.TRUE.)
              CALL PGCLOS(IDNEW)
              WRITE(*,101) '   ...OK! Plot created and closed'
              CALL PGSLCT(ID)
            END IF
            CALL BUTTON(NB,'postscript',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.23)THEN
            CALL BUTTON(NB,'RGB *.ppm',5)
            CALL CREATERGB(NCBUFF)
            CALL BUTTON(NB,'RGB *.ppm',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.24)THEN
            CALL BUTTON(NB,'[Cc]entroid',5)
            IF(CH.EQ.'C')THEN
              COFFSET(1:1)=READC('Measure relative offsets in box-9 (y/n)','n','yn')
            ELSE
              COFFSET='n'
            END IF
            IF(COFFSET.EQ.'n')THEN
              CALL CENTROID(NCBUFF,CH,XC_,YC_,.FALSE.,X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE, &
               EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,.TRUE.)
            ELSE
              NSIMUL_=NSIMUL
              NSIMUL=0
              LOK=.TRUE.
              LMIDE=.TRUE.
              DO WHILE(LMIDE)
                WRITE(*,100) 'Select object in frame #1...'
                CALL PGBAND(7,0,0.,0.,XC_,YC_,CH)
                IF(CH.NE.'X')THEN
                  WRITE(*,101) 'OK!'
                ELSE
                  WRITE(*,101) 'operation cancelled!'
                  LMIDE=.FALSE.
                  LOK=.FALSE.
                END IF
                IF(LOK)THEN
                  CALL CENTROID(NCBUFF,'d',XC_,YC_,.FALSE.,X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE, &
                   EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,.FALSE.)
                  WRITE(*,100) 'Please confirm fit with mouse...'
                  CALL PGBAND(0,0,0.,0.,XDUMMY,YDUMMY,CH)
                  IF(CH.NE.'X')THEN
                    LMIDE=.FALSE.
                    WRITE(*,101) 'OK!'
                  ELSE
                    WRITE(*,100) 'last fit has been rejected!'
                    WRITE(*,101) ' Try again.'
                  END IF
                END IF
              END DO
              IF(LOK)THEN
                XC_=X0-REAL(DJ(1)*NSIZEFB9(NCBUFF))
                YC_=Y0-REAL(DI(1)*NSIZEFB9(NCBUFF))
                XOFFSET(1)=0.
                YOFFSET(1)=0.
                DO I=2,NFRAMES(NCBUFF)
                  IF(LOK)THEN
                    LMIDE=.TRUE.
                    DO WHILE(LMIDE)
                      WRITE(*,'(A,I1,A,$)') 'Select the same object in frame #',I,'...'
                      CALL PGBAND(7,0,0.,0.,XC__,YC__,CH)
                      IF(CH.NE.'X')THEN
                        WRITE(*,101) 'OK!'
                      ELSE
                        WRITE(*,101) 'operation cancelled!'
                        LOK=.FALSE.
                        LMIDE=.FALSE.
                      END IF
                      IF(LMIDE)THEN
                        CALL CENTROID(NCBUFF,'d',XC__,YC__,.FALSE.,X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE, &
                         EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,.FALSE.)
                        WRITE(*,100) 'Please confirm fit with mouse...'
                        CALL PGBAND(0,0,0.,0.,XDUMMY,YDUMMY,CH)
                        IF(CH.NE.'X')THEN
                          LMIDE=.FALSE.
                          WRITE(*,101) 'OK!'
                        ELSE
                          WRITE(*,100) 'last fit has been rejected!'
                          WRITE(*,101) ' Try again.'
                        END IF
                      END IF
                    END DO
                    IF(LOK)THEN
                      XC__=X0-REAL(DJ(I)*NSIZEFB9(NCBUFF))
                      YC__=Y0-REAL(DI(I)*NSIZEFB9(NCBUFF))
                      XOFFSET(I)=XC__-XC_
                      YOFFSET(I)=YC__-YC_
                      WRITE(*,*) I,XOFFSET(I),YOFFSET(I)
                    END IF
                  END IF
                END DO
                IF(LOK)THEN
                  OFFFILE=READC('Output file name with offsets','@','@')
                  OPEN(10,FILE=OFFFILE,STATUS='UNKNOWN',FORM='FORMATTED')
                  DO I=1,NFRAMES(NCBUFF)
                    WRITE(10,*) I,XOFFSET(I),0.0,YOFFSET(I),0.0
                  END DO
                  CLOSE(10)
                END IF
              END IF
              NSIMUL=NSIMUL_
            END IF
            CALL BUTTON(NB,'[Cc]entroid',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.25)THEN
            CALL BUTTON(NB,'[r]ectify',5)
            IF(LPBOUND)THEN
              LPBOUND=.FALSE.
              CALL BUTTON(13,'plot bou.',0)
            END IF
            IF(LPMAPP)THEN
              LPMAPP=.FALSE.
              CALL BUTTON(18,'[p]lot mapp.',0)
            END IF
            NEWBUFF=NCBUFF
            DO WHILE((NEWBUFF.EQ.NCBUFF).AND.(NEWBUFF.NE.0))
              NEWBUFF=NCBUFF+1
              IF(NEWBUFF.GT.NMAXBUFF/2) NEWBUFF=1
              WRITE(CDUMMY,*) NEWBUFF
              NEWBUFF=READILIM('Buffer # to store rectified image (0=none)',CDUMMY,0,NMAXBUFF)
              IF(NEWBUFF.EQ.NCBUFF)THEN
                WRITE(*,101) '***ERROR***'
                WRITE(*,101) '=> New buffer must be different than current buffer. Try again.'
              END IF
            END DO
            IF(NEWBUFF.NE.0)THEN
              IF(IRESAMPL.EQ.1)THEN
                IRESAMPL=READILIM('Resampling mode (1=approx., 2=accurate)','1',1,2)
              ELSE
                IRESAMPL=READILIM('Resampling mode (1=approx., 2=accurate)','2',1,2)
              END IF
              CALL RECTIFY(NEWBUFF)
            END IF
            CALL BUTTON(NB,'[r]ectify',0)
            CALL BUTTON(NB,'[r]ectify',NCOLOR6)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.26)THEN
            CALL BUTTON(NB,'s[k]y sub.',5)
            IF(LPBOUND)THEN
              LPBOUND=.FALSE.
              CALL BUTTON(13,'plot bou.',0)
            END IF
            IF(LPMAPP)THEN
              LPMAPP=.FALSE.
              CALL BUTTON(18,'[p]lot mapp.',0)
            END IF
            LASK=.TRUE.
            DO WHILE(LASK)
              NEWBUFF=NCBUFF+1
              IF(NEWBUFF.GT.NMAXBUFF/2) NEWBUFF=1
              WRITE(CDUMMY,*) NEWBUFF
              NEWBUFF=READILIM('Buffer # to store intermediate data',CDUMMY,0,NMAXBUFF)
              LASK=(NEWBUFF.EQ.NCBUFF)
              IF(LASK)THEN
                WRITE(*,101) '***ERROR***'
                WRITE(*,101) '=> New buffer cannot be current buffer.'
                WRITE(*,100) '(press <CR> to continue...)'
                READ(*,*)
              END IF
            END DO
            LASK=.TRUE.
            DO WHILE(LASK)
              NEWBUFF_=NEWBUFF+1
              IF(NEWBUFF_.GT.NMAXBUFF/2) NEWBUFF_=1
              WRITE(CDUMMY,*) NEWBUFF_
              NEWBUFF_=READILIM('Buffer #2 to store intermediate data',CDUMMY,0,NMAXBUFF)
              LASK=(NEWBUFF_.EQ.NCBUFF)
              IF(LASK)THEN
                WRITE(*,101) '***ERROR***'
                WRITE(*,101) '=> New buffer cannot be current buffer.'
                WRITE(*,100) '(press <CR> to continue...)'
                READ(*,*)
              END IF
              LASK=(NEWBUFF_.EQ.NEWBUFF)
              IF(LASK)THEN
                WRITE(*,101) '***ERROR***'
                WRITE(*,101) '=> New buffer#2 cannot be new buffer #1.'
                WRITE(*,100) '(press <CR> to continue...)'
                READ(*,*)
              END IF
            END DO
            IF(IRESAMPL.EQ.1)THEN
              IRESAMPL=READILIM('Resampling mode (1=approx., 2=accurate)','1',1,2)
            ELSE
              IRESAMPL=READILIM('Resampling mode (1=approx., 2=accurate)','2',1,2)
            END IF
            CALL SKYSUB(NEWBUFF,NEWBUFF_)
            CALL BUTTON(NB,'s[k]y sub.',0)
            CALL BUTTON(NB,'s[k]y sub.',NCOLOR6)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.27)THEN
            CALL BUTTON(NB,'X-[e]xtract',5)
            CALL XEXTRACT
            CALL BUTTON(NB,'X-[e]xtract',0)
            CALL BUTTON(NB,'X-[e]xtract',NCOLOR6)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.28)THEN
            CALL BUTTON(NB,'[u]nrectify',5)
            I1=SYMINGRID-SYMINEXTG !borde inferior del primer pixel
            I2=SYMAXGRID+SYMAXEXTG !borde superior del ultimo pixel
            J1=SXMINGRID-SXMINEXTG !borde izquierdo del primer pixel
            J2=SXMAXGRID+SXMAXEXTG !borde derecho del ultimo pixel
            IF((NAXIS(1,NCBUFF).NE.(J2-J1)).OR.(NAXIS(2,NCBUFF).NE.(I2-I1)))THEN
              WRITE(*,101) '***ERROR***'
              WRITE(*,100) '=> image in current buffer does not have'
              WRITE(*,101) ' the expected dimensions.'
            ELSE
              LASK=.TRUE.
              DO WHILE(LASK)
                NEWBUFF=NCBUFF+1
                IF(NEWBUFF.GT.NMAXBUFF/2) NEWBUFF=1
                WRITE(CDUMMY,*) NEWBUFF
                NEWBUFF=READILIM('Buffer # to store intermediate image',CDUMMY,0,NMAXBUFF)
                LASK=(NEWBUFF.EQ.NCBUFF)
                IF(LASK)THEN
                  WRITE(*,101) '***ERROR***'
                  WRITE(*,101) '=> New buffer cannot be current buffer.'
                  WRITE(*,100) '(press <CR> to continue...)'
                  READ(*,*)
                END IF
              END DO
              IF(IRESAMPL.EQ.1)THEN
                IRESAMPL=READILIM('Resampling mode (1=approx., 2=accurate)','1',1,2)
              ELSE
                IRESAMPL=READILIM('Resampling mode (1=approx., 2=accurate)','2',1,2)
              END IF
              CALL UNRECTIFY(NEWBUFF)
            END IF
            CALL BUTTON(NB,'[u]nrectify',0)
            CALL BUTTON(NB,'[u]nrectify',NCOLOR6)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.29)THEN
            CALL BUTTON(NB,'[.]special',5)
            WRITE(*,*)
            WRITE(*,101) '(1) run SExtractor'
            WRITE(*,101) '(2) estimate slit aperture correction'
            WRITE(*,101) '(3) shift one image relative to another'
            WRITE(*,101) '(4) overplot ds9 regions'
            WRITE(*,101) '(5) overplot X,Y columns from ASCII file'
            WRITE(*,101) '(0) exit'
            CSPECIAL(1:1)=READC('Option (0..5)','0','012345')
            IF(CSPECIAL.EQ.'0')THEN
            ELSEIF(CSPECIAL.EQ.'1')THEN
              CALL SEXTRACTOR(NCBUFF)
            ELSEIF(CSPECIAL.EQ.'2')THEN
              CALL SLIT_APERTURE(NCBUFF)
            ELSEIF(CSPECIAL.EQ.'3')THEN
              CALL SHIFTIMAGE
            ELSEIF(CSPECIAL.EQ.'4')THEN
              LOGFILE=.FALSE.
              LOGFILERR=.FALSE.
              DO WHILE(.NOT.LOGFILE)
                INFILE_=READC('Input ds9 file name (none=EXIT, NONE=delete region)','*.reg','@')
                IF((INDEX(INFILE_,'*').NE.0).OR.(INDEX(INFILE_,'?').NE.0))THEN
                  L1=TRUEBEG(INFILE_)
                  L2=TRUELEN(INFILE_)
                  ISYSTEM=SYSTEMFUNCTION('ls '//INFILE_(L1:L2))
                ELSEIF(INFILE_.EQ.'none')THEN
                  LOGFILE=.TRUE.
                  LOGFILERR=.TRUE.
                ELSEIF(INFILE_.EQ.'NONE')THEN
                  LOGFILE=.TRUE.
                  LOGFILERR=.TRUE.
                  LDS9REG(NCBUFF)=.FALSE.
                ELSE
                  INQUIRE(FILE=INFILE_,EXIST=LOGFILE)
                  IF(.NOT.LOGFILE)THEN
                    !Try adding extension .reg
                    L1=TRUEBEG(INFILE_)
                    L2=TRUELEN(INFILE_)
                    INFILE_=INFILE_(L1:L2)//'.reg'
                    INQUIRE(FILE=INFILE_,EXIST=LOGFILE)
                    IF(.NOT.LOGFILE)THEN
                      WRITE(*,101) 'This file does not exist. Try again'
                      WRITE(*,101) INFILE_(L1:L2)
                    END IF
                  END IF
                END IF
              END DO
              IF((LOGFILE).AND.(.NOT.LOGFILERR))THEN
                LDS9REG(NCBUFF)=.TRUE.
                DS9REGFILE(NCBUFF)=INFILE_
                CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
              END IF
            ELSEIF(CSPECIAL.EQ.'5')THEN
              LOGFILE=.FALSE.
              LOGFILERR=.FALSE.
              DO WHILE(.NOT.LOGFILE)
                IF(ASCREGFILE(NCBUFF).NE.'none')THEN
                  CDUMMY=ASCREGFILE(NCBUFF)
                ELSE
                  CDUMMY='*'
                END IF
                INFILE_=READC('Input ASCII file (none=EXIT, NONE=delete region)',CDUMMY,'@')
                IF((INDEX(INFILE_,'*').NE.0).OR.(INDEX(INFILE_,'?').NE.0))THEN
                  L1=TRUEBEG(INFILE_)
                  L2=TRUELEN(INFILE_)
                  ISYSTEM=SYSTEMFUNCTION('ls '//INFILE_(L1:L2))
                ELSEIF(INFILE_.EQ.'none')THEN
                  LOGFILE=.TRUE.
                  LOGFILERR=.TRUE.
                ELSEIF(INFILE_.EQ.'NONE')THEN
                  LOGFILE=.TRUE.
                  LOGFILERR=.TRUE.
                  LASCREG(NCBUFF)=.FALSE.
                ELSE
                  INQUIRE(FILE=INFILE_,EXIST=LOGFILE)
                  IF(.NOT.LOGFILE)THEN
                    L1=TRUEBEG(INFILE_)
                    L2=TRUELEN(INFILE_)
                    WRITE(*,101) 'This file does not exist. Try again'
                    WRITE(*,101) INFILE_(L1:L2)
                  END IF
                END IF
              END DO
              IF((LOGFILE).AND.(.NOT.LOGFILERR))THEN
                WRITE(CDUMMY,*)NCOLASC1(NCBUFF)
                NCOLASC1(NCBUFF)=READI('Column number for X value',CDUMMY)
                WRITE(CDUMMY,*)NCOLASC2(NCBUFF)
                NCOLASC2(NCBUFF)=READI('Column number for Y value',CDUMMY)
                WRITE(CDUMMY,*)ASCSYMB(NCBUFF)
                ASCSYMB(NCBUFF)=READI('PGPLOT symbol number',CDUMMY)
                WRITE(CDUMMY,*)ASCCOLOR(NCBUFF)
                ASCCOLOR(NCBUFF)=READI('Color number',CDUMMY)
                WRITE(CDUMMY,*)ASCLWIDTH(NCBUFF)
                ASCLWIDTH(NCBUFF)=READI('Line width',CDUMMY)
                WRITE(CDUMMY,*)ASCNUMBER(NCBUFF)
                ASCNUMBER(NCBUFF)=READILIM('Plot number (0=no, 1=yes)',CDUMMY,0,1)
                WRITE(CDUMMY,*)ASCCHEIGHT(NCBUFF)
                ASCCHEIGHT(NCBUFF)=READF('Symbol size',CDUMMY)
                LASCREG(NCBUFF)=.TRUE.
                ASCREGFILE(NCBUFF)=INFILE_
                CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
              END IF
            END IF
            CALL BUTTON(NB,'[.]special',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.30)THEN
            CALL BUTTON(NB,'offsets',5)
            COFFSETS(1:1)=READC('Use [c]entroid or just c[u]rsor',COFFSETS,'cu')
            LOGFILE=.FALSE.
            DO WHILE(.NOT.LOGFILE)
              FILELISTIN=READC('Input file name with list of FITS files','@','@')
              INQUIRE(FILE=FILELISTIN,EXIST=LOGFILE)
              IF(.NOT.LOGFILE)THEN
                WRITE(*,100) 'ERROR: this file does not exist.'
                WRITE(*,101) ' Try again.'
              END IF
            END DO
            LOGFILE=.TRUE.
            DO WHILE(LOGFILE)
              FILELISTOUT=READC('Output file name with list of FITS files','@','@')
              INQUIRE(FILE=FILELISTOUT,EXIST=LOGFILE)
              IF(LOGFILE)THEN
                WRITE(*,100) 'ERROR: this file already exists.'
                WRITE(*,101) ' Try again.'
              END IF
            END DO
            OPEN(10,FILE=FILELISTIN,STATUS='OLD',FORM='FORMATTED')
            NFITSFILES=0
            L3=0
880         READ(10,101,END=882) FILEOFFSET
            INQUIRE(FILE=FILEOFFSET,EXIST=LOGFILE)
            IF(.NOT.LOGFILE)THEN
              WRITE(*,101) 'ERROR: the following file does not exist:'
              WRITE(*,101) FILEOFFSET(1:TRUELEN(FILEOFFSET))
              GOTO 888
            END IF
            NFITSFILES=NFITSFILES+1
            L1=TRUEBEG(FILEOFFSET)
            L2=TRUELEN(FILEOFFSET)
            IF(L2-L1+1.GT.L3) L3=L2-L1+1
            GOTO 880
882         CLOSE(10)
            IF(NFITSFILES.LE.1)THEN
              WRITE(*,101) 'ERROR: insufficient number of fits files.'
              GOTO 888
            END IF
            OPEN(10,FILE=FILELISTIN,STATUS='OLD',FORM='FORMATTED')
            OPEN(12,FILE=FILELISTOUT,STATUS='NEW',FORM='FORMATTED')
            DO I=1,NFITSFILES
              READ(10,101) FILEOFFSET
              L1=TRUEBEG(FILEOFFSET)
              L2=TRUELEN(FILEOFFSET)
              CALL SLEEFITS(FILEOFFSET(L1:L2),.FALSE.,IROTATE,NCBUFF,.FALSE.,0)
              NFRAMES(NCBUFF)=1
              NSIZEFB9(NCBUFF)=0
              NAXISFRAME(1,1,NCBUFF)=0
              NAXISFRAME(2,1,NCBUFF)=0
              INFILE(NCBUFF)=FILEOFFSET
              DO K=1,NMAXBUFF
                IF(K.LE.NMAXBUFF/2)THEN
                  WRITE(CDUMMY,'(A2,I1,A1)') '#[',K,']'
                  CALL RMBLANK(CDUMMY,CDUMMY,L)
                  IF(K.EQ.NCBUFF)THEN
                    CALL BUTTON(34+K,CDUMMY(1:L),1)
                  ELSE
                    CALL BUTTON(34+K,CDUMMY(1:L),0)
                  END IF
                ELSE
                  WRITE(CDUMMY,'(A5,I1)') 'err #',K-NMAXBUFF/2
                  L=TRUELEN(CDUMMY)
                  IF(K.EQ.NCBUFF)THEN
                    CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),1)
                  ELSE
                    CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),0)
                  END IF
                END IF
              END DO
              NX1=1
              NX2=NAXIS(1,NCBUFF)
              NY1=1
              NY2=NAXIS(2,NCBUFF)
              CALL STATISTIC(NCBUFF,NX1,NX2,NY1,NY2,.FALSE.,.TRUE.,.TRUE.,0.0,.FALSE.)
              IF(LFIRSTPLOT)THEN
                CALL ACTIVEBUT(JUST,MODECUT)
                IF(FSIGMA.GT.0.0)THEN
                  BG=FMEAN-5.*FSIGMA
                  FG=FMEAN+5.*FSIGMA
                ELSE
                  BG=FMEAN-1.0
                  FG=FMEAN+1.0
                END IF
                CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
                LFIRSTPLOT=.FALSE.
              ELSE
                CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
              END IF
              LREPEAT=.TRUE.
              DO WHILE(LREPEAT)
                WRITE(*,100) 'Select object in image'
                WRITE(*,100) FILEOFFSET(L1:L2)
                WRITE(*,100) '...'
                CALL PGBAND(7,0,0.,0.,XC_,YC_,CH)
                IF(CH.NE.'X')THEN
                  WRITE(*,101) 'OK!'
                ELSE
                  WRITE(*,101) 'operation cancelled!'
                  CLOSE(10)
                  CLOSE(12)
                  GOTO 888
                END IF
                IF(COFFSETS.EQ.'c')THEN
                  CALL CENTROID(NCBUFF,'d',XC_,YC_,.FALSE.,X0,Y0,SIGMAX,SIGMAY,BETA,AMP,CTE, &
                   EX0,EY0,ESIGMAX,ESIGMAY,EBETA,EAMP,ECTE,.FALSE.)
                  WRITE(*,100) 'Please confirm fit with mouse...'
                ELSE
                  X0 = XC_
                  Y0 = YC_
                  WRITE(*,100) 'Please confirm with mouse...'
                END IF
                CALL PGBAND(0,0,0.,0.,XDUMMY,YDUMMY,CH)
                IF(CH.NE.'X')THEN
                  LREPEAT=.FALSE.
                  WRITE(*,101) 'OK!'
                ELSE
                  IF(COFFSETS.EQ.'c')THEN
                    WRITE(*,100) 'last fit has been rejected!'
                  ELSE
                    WRITE(*,100) 'last position has been rejected!'
                  END IF
                  WRITE(*,101) ' Try again.'
                END IF
              END DO
              IF(I.EQ.1)THEN
                XOFFSET_ORIGEN=X0
                YOFFSET_ORIGEN=Y0
                X0=0
                Y0=0
              ELSE
                X0=XOFFSET_ORIGEN-X0
                Y0=YOFFSET_ORIGEN-Y0
              END IF
              WRITE(12,100) '0  '  ! include flag for imcombine program
              WRITE(12,100) FILEOFFSET(L1:L2)
              IF(L2-L1+1.LT.L3)THEN
                DO L=1,L3-(L2-L1+1)
                  WRITE(12,100) ' '
                END DO
              END IF
              WRITE(12,*) NINT(X0),NINT(Y0)
            END DO
            CLOSE(10)
            CLOSE(12)
888         CALL BUTTON(NB,'offsets',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.31)THEN
            CALL BUTTON(NB,'resize ima',5)
            WRITE(*,100) 'NAXIS(1) (x-axis): '
            WRITE(*,*) NAXIS(1,NCBUFF)
            WRITE(*,100) 'NAXIS(2) (y-axis): '
            WRITE(*,*) NAXIS(2,NCBUFF)
            WRITE(CDUMMY,*) NAXIS(1,NCBUFF)
            NAXIS1_=READILIM('New NAXIS(1)',CDUMMY,1,NXMAX)
            WRITE(CDUMMY,*) NAXIS(2,NCBUFF)
            NAXIS2_=READILIM('New NAXIS(2)',CDUMMY,1,NYMAX)
            DAXIS1_=READI('Integer offset for X shift (+right,-left)','0')
            DAXIS2_=READI('Integer offset for Y shift ...(+up,-down)','0')
! inicializamos imagen de destino a cero
            DO I=1,NAXIS2_
              DO J=1,NAXIS1_
                IMAGEN_(J,I)=0.0
              END DO
            END DO
! trasladamos los datos a la imagen de destino
            DO I=1,NAXIS(2,NCBUFF)
              II=I+DAXIS2_
              IF((II.GE.1).AND.(II.LE.NAXIS2_))THEN
                DO J=1,NAXIS(1,NCBUFF)
                  JJ=J+DAXIS1_
                  IF((JJ.GE.1).AND.(JJ.LE.NAXIS1_))THEN
                    IMAGEN_(JJ,II)=IMAGEN(J,I,NCBUFF)
                  END IF
                END DO
              END IF
            END DO
! pasamos los datos de la matriz temporal a la definitiva
            DO I=1,NAXIS2_
              DO J=1,NAXIS1_
                IMAGEN(J,I,NCBUFF)=IMAGEN_(J,I)
              END DO
            END DO
! Repetimos el trabajo con la imagen de errores
            IF(NCBUFF.LE.NMAXBUFF/2)THEN
              DO I=1,NAXIS2_
                DO J=1,NAXIS1_
                  IMAGEN_(J,I)=0.0
                END DO
              END DO
              DO I=1,NAXIS(2,NCBUFF)
                II=I+DAXIS2_
                IF((II.GE.1).AND.(II.LE.NAXIS2_))THEN
                  DO J=1,NAXIS(1,NCBUFF)
                    JJ=J+DAXIS1_
                    IF((JJ.GE.1).AND.(JJ.LE.NAXIS1_))THEN
                      IMAGEN_(JJ,II)=IMAGEN(J,I,NCBUFF+NMAXBUFF/2)
                    END IF
                  END DO
                END IF
              END DO
              DO I=1,NAXIS2_
                DO J=1,NAXIS1_
                  IMAGEN(J,I,NCBUFF+NMAXBUFF/2)=IMAGEN_(J,I)
                END DO
              END DO
            END IF
! redefinimos las dimensiones
            NAXIS(1,NCBUFF)=NAXIS1_
            NAXIS(2,NCBUFF)=NAXIS2_
            IF(NCBUFF.LE.NMAXBUFF/2)THEN
              NAXIS(1,NCBUFF+NMAXBUFF/2)=NAXIS1_
              NAXIS(2,NCBUFF+NMAXBUFF/2)=NAXIS2_
            END IF
            WRITE(*,101) 'WARNING: replot image to see updated size.'
            CALL BUTTON(NB,'resize ima',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.32)THEN
            CALL BUTTON(NB,'combine',5)
            IF(NFRAMES(NCBUFF).GT.1)THEN
              BG_=BG
              FG_=FG
              NX1_=NX1
              NX2_=NX2
              NY1_=NY1
              NY2_=NY2
              CALL COMBINE(NCBUFF)
              BG=BG_
              FG=FG_
              NX1=NX1_
              NX2=NX2_
              NY1=NY1_
              NY2=NY2_
              CALL HISTOGRAM(NCBUFF)
              CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              IF(LSTACK.AND.(NSTACK.GT.0))THEN
                CALL PGSCI(3)
                CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
                CALL PGSCI(1)
              END IF
              IF(LPBOUND) CALL PBOUND
              IF(LPMAPP) CALL PMAPP
              IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            ELSE
              WRITE(*,100) '> NFRAMES(NCBUFF): '
              WRITE(*,*) NFRAMES(NCBUFF)
              WRITE(*,101) 'ERROR: insufficient number of frames/image'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
            END IF
            CALL BUTTON(NB,'combine',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.33)THEN
            CALL BUTTON(NB,'shiftb9',5)
            IF(NFRAMES(NCBUFF).GT.1)THEN
              BG_=BG
              FG_=FG
              NX1_=NX1
              NX2_=NX2
              NY1_=NY1
              NY2_=NY2
              CALL SHIFTB9(NCBUFF)
              BG=BG_
              FG=FG_
              NX1=NX1_
              NX2=NX2_
              NY1=NY1_
              NY2=NY2_
              CALL HISTOGRAM(NCBUFF)
              CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              IF(LSTACK.AND.(NSTACK.GT.0))THEN
                CALL PGSCI(3)
                CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
                CALL PGSCI(1)
              END IF
              IF(LPBOUND) CALL PBOUND
              IF(LPMAPP) CALL PMAPP
              IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            ELSE
              WRITE(*,100) '> NFRAMES(NCBUFF): '
              WRITE(*,*) NFRAMES(NCBUFF)
              WRITE(*,101) 'ERROR: insufficient number of frames/image'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
            END IF
            CALL BUTTON(NB,'shiftb9',0)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.34)THEN
            CALL BUTTON(NB,'box[9]oper',5)
            IF(NBUFFBOX9.EQ.0)THEN
              NBUFFBOX9=READILIM('Buffer with set of box-9 frames (0=exit)','@',0,NMAXBUFF)
            ELSE
              WRITE(CDUMMY,*) NBUFFBOX9
              NBUFFBOX9=READILIM('Buffer with set of box-9 frames (0=exit)',CDUMMY,0,NMAXBUFF)
            END IF
            IF(NBUFFBOX9.GT.0)THEN
              BG_=BG
              FG_=FG
              NX1_=NX1
              NX2_=NX2
              NY1_=NY1
              NY2_=NY2
              CALL OPERBOX9(NBUFFBOX9)
              BG=BG_
              FG=FG_
              NX1=NX1_
              NX2=NX2_
              NY1=NY1_
              NY2=NY2_
              CALL HISTOGRAM(NCBUFF)
              CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              IF(LSTACK.AND.(NSTACK.GT.0))THEN
                CALL PGSCI(3)
                CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
                CALL PGSCI(1)
              END IF
              IF(LPBOUND) CALL PBOUND
              IF(LPMAPP) CALL PMAPP
              IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            END IF
            CALL BUTTON(NB,'box[9]oper',0)
!------------------------------------------------------------------------------
          ELSEIF((NB.GE.35).AND.(NB.LE.34+NMAXBUFF/2))THEN
            IF((NAXIS(1,NB-34).EQ.0).AND.(NAXIS(2,NB-34).EQ.0))THEN
              WRITE(*,*)
              WRITE(*,101) '***WARNING***'
              WRITE(*,100) '=> Buffer #'
              WRITE(*,*) NB-34
              WRITE(*,101) '=> This DATA buffer has not been defined.'
            ELSE
              NCBUFF_=NCBUFF
              NCBUFF=NB-34
              DO K=1,NMAXBUFF/2
                WRITE(CDUMMY,'(A2,I1,A1)') '#[',K,']'
                CALL RMBLANK(CDUMMY,CDUMMY,L)
                IF(K.EQ.NCBUFF)THEN
                  CALL BUTTON(34+K,CDUMMY(1:L),1)
                ELSE
                  CALL BUTTON(34+K,CDUMMY(1:L),0)
                END IF
              END DO
              DO K=NMAXBUFF/2+1,NMAXBUFF
                WRITE(CDUMMY,'(A5,I1)') 'err #',K-NMAXBUFF/2
                L=TRUELEN(CDUMMY)
                CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),0)
              END DO
              IF((NAXIS(1,NCBUFF).EQ.NAXIS(1,NCBUFF_)).AND.(NAXIS(2,NCBUFF).EQ.NAXIS(2,NCBUFF_)))THEN
                CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
              ELSE
                NX1=1
                NX2=NAXIS(1,NCBUFF)
                NY1=1
                NY2=NAXIS(2,NCBUFF)
                CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              END IF
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.41)THEN
            CALL BUTTON(NB,'zoom',5)
            CALL HPLOT('zoom')
            CALL BUTTON(NB,'zoom',0)
            CALL BUTTON(NB,'zoom',NCOLOR4)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.42)THEN
            CALL BUTTON(NB,'whole',5)
            CALL HPLOT('whole')
            CALL BUTTON(NB,'whole',0)
            CALL BUTTON(NB,'whole',NCOLOR4)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.43)THEN
            CALL BUTTON(NB,'min,max',5)
            CALL HPLOT('min,max')
            CALL BUTTON(NB,'min,max',0)
            CALL BUTTON(NB,'min,max',NCOLOR4)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.44)THEN
            CALL BUTTON(NB,'[o]ffsets',5)
            IF(NBUFFBOX9.EQ.0)THEN
              NBUFFBOX9=READILIM('Buffer with frames to be measured (0=exit)','@',0,NMAXBUFF)
            ELSE
              WRITE(CDUMMY,*) NBUFFBOX9
              NBUFFBOX9=READILIM('Buffer with frames to be measured (0=exit)',CDUMMY,0,NMAXBUFF)
            END IF
            IF(NBUFFBOX9.GT.0)THEN
              BG_=BG
              FG_=FG
              NX1_=NX1
              NX2_=NX2
              NY1_=NY1
              NY2_=NY2
              CALL MIDEOFF(NBUFFBOX9)
              BG=BG_
              FG=FG_
              NX1=NX1_
              NX2=NX2_
              NY1=NY1_
              NY2=NY2_
              CALL HISTOGRAM(NCBUFF)
              CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              IF(LSTACK.AND.(NSTACK.GT.0))THEN
                CALL PGSCI(3)
                CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
                CALL PGSCI(1)
              END IF
              IF(LPBOUND) CALL PBOUND
              IF(LPMAPP) CALL PMAPP
              IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            END IF
            CALL BUTTON(NB,'[o]ffsets',0)
!------------------------------------------------------------------------------
          ELSEIF((NB.GE.45).AND.(NB.LE.44+NMAXBUFF/2))THEN
            IF((NAXIS(1,NB-44+NMAXBUFF/2).EQ.0).AND.(NAXIS(2,NB-44+NMAXBUFF/2).EQ.0))THEN
              WRITE(*,101) '***WARNING***'
              WRITE(*,101) '=> This ERROR buffer has not been defined.'
            ELSE
              NCBUFF_=NCBUFF
              NCBUFF=NB-44+NMAXBUFF/2
              DO K=NMAXBUFF/2+1,NMAXBUFF
                WRITE(CDUMMY,'(A5,I1)') 'err #',K-NMAXBUFF/2
                L=TRUELEN(CDUMMY)
                IF(K.EQ.NCBUFF)THEN
                  CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),1)
                ELSE
                  CALL BUTTON(44+K-NMAXBUFF/2,CDUMMY(1:L),0)
                END IF
              END DO
              DO K=1,NMAXBUFF/2
                WRITE(CDUMMY,'(A2,I1,A1)') '#[',K,']'
                L=TRUELEN(CDUMMY)
                CALL BUTTON(34+K,CDUMMY(1:L),0)
              END DO
              IF((NAXIS(1,NCBUFF).EQ.NAXIS(1,NCBUFF_)).AND.(NAXIS(2,NCBUFF).EQ.NAXIS(2,NCBUFF_)))THEN
                CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
              ELSE
                NX1=1
                NX2=NAXIS(1,NCBUFF)
                NY1=1
                NY2=NAXIS(2,NCBUFF)
                CALL SUBLOOK(.FALSE.,NCBUFF,.FALSE.)
              END IF
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.51)THEN
            IF(LOVERCUTS)THEN
              CALL BUTTON(51,'over=FALSE',0)
              CALL BUTTON(51,'over=FALSE',NCOLOR4)
              LOVERCUTS=.FALSE.
            ELSE
              CALL BUTTON(51,'over=TRUE',0)
              CALL BUTTON(51,'over=TRUE',NCOLOR4)
              LOVERCUTS=.TRUE.
            END IF
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.161)THEN
            CALL BUTTON(NB,'zoom',5)
            CALL ZOOMHISTOG(BG,FG)
            CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
            IF(LSTACK.AND.(NSTACK.GT.0))THEN
              CALL PGSCI(3)
              CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
              CALL PGSCI(1)
            END IF
            IF(LPBOUND) CALL PBOUND
            IF(LPMAPP) CALL PMAPP
            IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            CALL STATISTIC(NCBUFF,NX1,NX2,NY1,NY2,.FALSE.,.FALSE.,.FALSE.,0.0,.FALSE.)
            CALL HISTOGRAM(NCBUFF)
            CALL BUTTON(NB,'zoom',0)
            CALL BUTTON(NB,'zoom',NCOLOR3)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.162)THEN
            CALL BUTTON(NB,'min[,]max',5)
            CALL STATISTIC(NCBUFF,NX1,NX2,NY1,NY2,.FALSE.,.FALSE.,.TRUE.,0.0,.FALSE.)
            IF(FMIN.NE.FMAX)THEN
              BG=FMIN
              FG=FMAX
            ELSE
              BG=FMIN-1.0
              FG=FMIN+1.0
            END IF
            CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
            IF(LSTACK.AND.(NSTACK.GT.0))THEN
              CALL PGSCI(3)
              CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
              CALL PGSCI(1)
            END IF
            IF(LPBOUND) CALL PBOUND
            IF(LPMAPP) CALL PMAPP
            IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            CALL HISTOGRAM(NCBUFF)
            CALL BUTTON(NB,'min[,]max',0)
            CALL BUTTON(NB,'min[,]max',NCOLOR3)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.163)THEN
            CALL BUTTON(NB,'z1[/]z2',5)
            IF(NFRAMES(NCBUFF).GT.1)THEN
              CALL STATISTICB9(NFRAMES(NCBUFF),NCBUFF,FMEAN,FSIGMA)
              IF(FSIGMA.GT.0.0)THEN
                BG=FMEAN-5.*FSIGMA
                FG=FMEAN+5.*FSIGMA
                WRITE(*,100) '>>> Mean : '
                WRITE(*,*) FMEAN
                WRITE(*,100) '>>> Sigma: '
                WRITE(*,*) FSIGMA
              ELSE
                BG=FMEAN-1.0
                FG=FMEAN+1.0
              END IF
            ELSE
              CALL ZSCALE(NCBUFF,NX1,NX2,NY1,NY2,Z1,Z2,FMEAN,FSIGMA)
              BG=Z1
              FG=Z2
            END IF
            CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
            IF(LSTACK.AND.(NSTACK.GT.0))THEN
              CALL PGSCI(3)
              CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
              CALL PGSCI(1)
            END IF
            IF(LPBOUND) CALL PBOUND
            IF(LPMAPP) CALL PMAPP
            IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            CALL HISTOGRAM(NCBUFF)
            CALL BUTTON(NB,'z1[/]z2',0)
            CALL BUTTON(NB,'z1[/]z2',NCOLOR3)
!------------------------------------------------------------------------------
          ELSEIF(NB.EQ.164)THEN
            CALL BUTTON(NB,'BG[:]FG',5)
            WRITE(CDUMMY,*) BG
            BG=READF('New BG',CDUMMY)
            WRITE(CDUMMY,*) FG
            FG=READF('New FG',CDUMMY)
            IF(FG.EQ.BG)THEN
              BG=BG-0.5
              FG=FG+0.5
            END IF
            CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
            IF(LSTACK.AND.(NSTACK.GT.0))THEN
              CALL PGSCI(3)
              CALL PGPOINT(NSTACK,XSTACK,YSTACK,21)
              CALL PGSCI(1)
            END IF
            IF(LPBOUND) CALL PBOUND
            IF(LPMAPP) CALL PMAPP
            IF(LPWAVE) CALL PWAVE(NDEGWV,CWV)
            CALL HISTOGRAM(NCBUFF)
            CALL BUTTON(NB,'BG[:]FG',0)
            CALL BUTTON(NB,'BG[:]FG',NCOLOR3)
!------------------------------------------------------------------------------
          END IF
!------------------------------------------------------------------------------
        END DO
!------------------------------------------------------------------------------
! ***END OF MAIN LOOP***
!------------------------------------------------------------------------------
! Close graphic display and end program
        CALL PGEND
        CALL Deallocate_Array_IMAGEN
        CALL Deallocate_Array_IMAGEN_
        STOP
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END PROGRAM XNIRSPEC
