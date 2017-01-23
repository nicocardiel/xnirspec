        SUBROUTINE SEXTRACTOR(NCBUFF)
        IMPLICIT NONE
        INTEGER NCBUFF
C
        INCLUDE 'dimensions.inc'
        INCLUDE 'sexpath.inc'
C
        INTEGER MAX_NCOMMENTS_CAT
        PARAMETER (MAX_NCOMMENTS_CAT=200)
C
        INTEGER TRUEBEG,TRUELEN
        INTEGER SYSTEMFUNCTION
        INTEGER READI
        INTEGER PGOPEN
        REAL READF
        REAL LEECOLUMN
        CHARACTER*255 READC
C
        INTEGER I,II,ISUM,N0,N0_
        INTEGER L1,L2,LL1,LL2,LLL1,LLL2
        INTEGER N1_,N2_
        INTEGER ISYSTEM
        INTEGER NCOMMENTS_CAT
        INTEGER N_X_IMAGE,N_Y_IMAGE
        INTEGER N_XMIN_IMAGE,N_XMAX_IMAGE,N_YMIN_IMAGE,N_YMAX_IMAGE
        INTEGER N_A_IMAGE,N_B_IMAGE,N_THETA_IMAGE
        INTEGER N_CLASS_STAR
        INTEGER NDATA_REGISTER(MAX_NCOMMENTS_CAT)
        INTEGER IDNEW,IDOLD
        REAL X_IMAGE,Y_IMAGE
        REAL XMIN_IMAGE,XMAX_IMAGE,YMIN_IMAGE,YMAX_IMAGE
        REAL A_IMAGE,B_IMAGE,THETA_IMAGE
        REAL CLASS_STAR
        REAL X0,Y0,R0,THETA0
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XELL(361),YELL(361)
        REAL FDUMMY
        REAL PHOTFLAM(NMAXBUFF),PHOTZPT(NMAXBUFF),EXPTIME(NMAXBUFF)
        REAL ZEROPOINT(NMAXBUFF)
        REAL CHOLD
        CHARACTER*1 COPC,CREF
        CHARACTER*22 CLABEL(MAX_NCOMMENTS_CAT)
        CHARACTER*30 FITSFILE,CMASK
        CHARACTER*50 CDUMMY
        CHARACTER*255 REFFILE,SEXPATH_
        CHARACTER*1000 CLINEA
        LOGICAL L_X_IMAGE,L_Y_IMAGE
        LOGICAL L_XMIN_IMAGE,L_XMAX_IMAGE
        LOGICAL L_YMIN_IMAGE,L_YMAX_IMAGE,L_A_IMAGE,L_B_IMAGE
        LOGICAL L_THETA_IMAGE,L_CLASS_STAR
        LOGICAL LOGFILE,LOGCATFILE
        LOGICAL LOOP,LOK
        LOGICAL L_PHOTFLAM(NMAXBUFF),L_PHOTZPT(NMAXBUFF)
        LOGICAL L_EXPTIME(NMAXBUFF)
        LOGICAL LPOSTSCRIPT
C
        COMMON/BLK_HSTFLUX1/L_PHOTFLAM,L_PHOTZPT,L_EXPTIME
        COMMON/BLK_HSTFLUX2/PHOTFLAM,PHOTZPT,EXPTIME
        COMMON/BLK_HSTFLUX3/ZEROPOINT
C------------------------------------------------------------------------------
        LLL1=0 !evita WARNING de compilacion
        LLL2=0 !evita WARNING de compilacion
        L_X_IMAGE=.FALSE. !evita WARNING de compilacion
        L_Y_IMAGE=.FALSE. !evita WARNING de compilacion
        L_XMIN_IMAGE=.FALSE. !evita WARNING de compilacion
        L_XMAX_IMAGE=.FALSE. !evita WARNING de compilacion
        L_YMIN_IMAGE=.FALSE. !evita WARNING de compilacion
        L_YMAX_IMAGE=.FALSE. !evita WARNING de compilacion
        L_A_IMAGE=.FALSE. !evita WARNING de compilacion
        L_B_IMAGE=.FALSE. !evita WARNING de compilacion
        L_THETA_IMAGE=.FALSE. !evita WARNING de compilacion
C
        SEXPATH_=SEXPATH
        CALL PGQWIN(XMIN,XMAX,YMIN,YMAX)
C Main loop
        LOOP=.TRUE.
        DO WHILE(LOOP)
          WRITE(*,*)
          INQUIRE(FILE='test.cat',EXIST=LOGCATFILE)
          IF(LOGCATFILE)THEN !vemos que cosas hay en el fichero
            L_X_IMAGE=.FALSE.
            L_Y_IMAGE=.FALSE.
            L_XMIN_IMAGE=.FALSE.
            L_XMAX_IMAGE=.FALSE.
            L_YMIN_IMAGE=.FALSE.
            L_YMAX_IMAGE=.FALSE.
            L_A_IMAGE=.FALSE.
            L_B_IMAGE=.FALSE.
            L_THETA_IMAGE=.FALSE.
            L_CLASS_STAR=.FALSE.
            OPEN(44,FILE='test.cat',STATUS='OLD',FORM='FORMATTED')
6           READ(44,101,END=8) CLINEA
            IF(CLINEA(1:1).EQ.'#')THEN
              IF(CLINEA(7:22).EQ.'X_IMAGE') L_X_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'Y_IMAGE') L_Y_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'XMIN_IMAGE') L_XMIN_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'XMAX_IMAGE') L_XMAX_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'YMIN_IMAGE') L_YMIN_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'YMAX_IMAGE') L_YMAX_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'A_IMAGE') L_A_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'B_IMAGE') L_B_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'THETA_IMAGE') L_THETA_IMAGE=.TRUE.
              IF(CLINEA(7:22).EQ.'CLASS_STAR') L_CLASS_STAR=.TRUE.
            ELSE
              GOTO 8
            END IF
            GOTO 6
8           CLOSE(44)
          END IF
          WRITE(*,101) '(R) Run SExtractor over current frame'
          IF(LOGCATFILE)THEN
            CMASK='Ri0x'
            WRITE(*,101) '(i) Info in a given object'
            WRITE(*,101) '(0) Plot number in catalogue'
            IF(L_X_IMAGE.AND.L_Y_IMAGE)THEN
              WRITE(*,101) '(1) Plot barycenter'
              CMASK=CMASK(1:TRUELEN(CMASK))//'1'
            END IF
            IF(L_XMIN_IMAGE.AND.L_XMAX_IMAGE.AND.
     +         L_YMIN_IMAGE.AND.L_YMAX_IMAGE)THEN
              WRITE(*,101) '(2) Plot rectangle which encloses '//
     +         'the detected object'
              CMASK=CMASK(1:TRUELEN(CMASK))//'2'
            END IF
            IF(L_A_IMAGE.AND.L_B_IMAGE.AND.L_THETA_IMAGE)THEN
              WRITE(*,101) '(3) Plot isophotal limit (R=3)'
              CMASK=CMASK(1:TRUELEN(CMASK))//'3'
            END IF
          END IF
          WRITE(*,101) '(x) EXIT'
          IF(LOGCATFILE)THEN
            COPC(1:1)=READC('Option (R/i/0..3/x)','x',CMASK)
          ELSE
            COPC(1:1)=READC('Option (R/x)','x','Rx')
          END IF
C
        IF(COPC.EQ.'0')THEN
          LPOSTSCRIPT=
     +     (READC('Create postscript file (y/n)','n','yn').EQ.'y')
        ELSE
          LPOSTSCRIPT=.FALSE.
        END IF
C..............................................................................
          IF(COPC.EQ.'x')THEN
            RETURN
C..............................................................................
          ELSEIF(COPC.EQ.'R')THEN
            LOK=.TRUE.
            INQUIRE(FILE='default.conv',EXIST=LOGFILE)
            IF(.NOT.LOGFILE)THEN
              WRITE(*,101) 'ERROR: file default.conv does not exist.'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
              LOK=.FALSE.
            ELSE
              WRITE(*,101) '> File default.conv exists. OK!'
            END IF
            INQUIRE(FILE='default.nnw',EXIST=LOGFILE)
            IF(.NOT.LOGFILE)THEN
              WRITE(*,101) 'ERROR: file default.nnw does not exist.'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
              LOK=.FALSE.
            ELSE
              WRITE(*,101) '> File default.nnw exists. OK!'
            END IF
            INQUIRE(FILE='default.param',EXIST=LOGFILE)
            IF(.NOT.LOGFILE)THEN
              WRITE(*,101) 'ERROR: file default.param does not exist.'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
              LOK=.FALSE.
            ELSE
              WRITE(*,101) '> File default.param exists. OK!'
            END IF
            INQUIRE(FILE='default.sex',EXIST=LOGFILE)
            IF(.NOT.LOGFILE)THEN
              WRITE(*,101) 'ERROR: file default.conv does not exist.'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
              LOK=.FALSE.
            ELSE
              WRITE(*,101) '> File default.sex exists. OK!'
            END IF
            CREF(1:1)=
     +       READC('Are you running SExtractor in double-image '//
     +       'mode (y/n)','n','yn')
            IF(CREF.EQ.'y')THEN
              LOGFILE=.FALSE.
              DO WHILE(.NOT.LOGFILE)
                REFFILE=READC('Reference image to be used for '//
     +           'detection of sources','*.fits','@')
                IF((INDEX(REFFILE,'*').NE.0).OR.
     +             (INDEX(REFFILE,'?').NE.0))THEN
                  ISYSTEM=SYSTEMFUNCTION('ls '//
     +             REFFILE(1:TRUELEN(REFFILE)))
                ELSE
                  INQUIRE(FILE=REFFILE,EXIST=LOGFILE)
                  IF(.NOT.LOGFILE)THEN
                    WRITE(*,101) 'ERROR: this file does not exist. '//
     +               'Try again.'
                    WRITE(*,100) 'Press <CR> to continue...'
                    READ(*,*)
                  ELSE
                    LLL1=TRUEBEG(REFFILE)
                    LLL2=TRUELEN(REFFILE)
                  END IF
                END IF
              END DO
            END IF
            IF(LOK)THEN
              INQUIRE(FILE='xyzxyzxyz.fits',EXIST=LOGFILE)
              IF(LOGFILE)THEN
                WRITE(*,101) '> Removing previous file xyzxyzxyz.fits'
                ISYSTEM=SYSTEMFUNCTION('rm -f xyzxyzxyz.fits')
              END IF
              WRITE(*,101) '> Creating new file xyzxyzxyz.fits'
              FITSFILE='xyzxyzxyz.fits'
              CALL SESCRFITS(FITSFILE,NCBUFF)
              L1=TRUEBEG(SEXPATH)
              L2=TRUELEN(SEXPATH)
              IF(L_PHOTFLAM(NCBUFF).AND.L_PHOTZPT(NCBUFF).AND.
     +         L_EXPTIME(NCBUFF))THEN  !HST flux calibration keywords available
                WRITE(CDUMMY,*) PHOTFLAM(NCBUFF)
                PHOTFLAM(NCBUFF)=READF('PHOTFLAM',CDUMMY)
                WRITE(CDUMMY,*) PHOTZPT(NCBUFF)
                PHOTZPT(NCBUFF)=READF('PHOTZPT',CDUMMY)
                WRITE(CDUMMY,*) EXPTIME(NCBUFF)
                EXPTIME(NCBUFF)=READF('EXPTIME',CDUMMY)
                ZEROPOINT(NCBUFF)=PHOTZPT(NCBUFF)
     +           -2.5*ALOG10(PHOTFLAM(NCBUFF)/EXPTIME(NCBUFF))
                WRITE(CDUMMY,*) ZEROPOINT(NCBUFF)
                ZEROPOINT(NCBUFF)=READF('MAG_ZEROPOINT',CDUMMY)
              ELSE
                ZEROPOINT(NCBUFF)=READF('MAG_ZEROPOINT','0.0')
              END IF
              WRITE(*,101) '> Running SExtractor...'
              IF(ZEROPOINT(NCBUFF).NE.0.0)THEN
                WRITE(CDUMMY,*) ZEROPOINT(NCBUFF)
                LL1=TRUEBEG(CDUMMY)
                LL2=TRUELEN(CDUMMY)
                IF(CREF.EQ.'y')THEN
                  WRITE(*,101) '% '//SEXPATH_(L1:L2)//'/sex \\'
                  WRITE(*,100) '  '//REFFILE(LLL1:LLL2)
                  WRITE(*,100) ',xyzxyzxyz.fits -MAG_ZEROPOINT '
                  WRITE(*,101) CDUMMY(LL1:LL2)
                  ISYSTEM=SYSTEMFUNCTION(SEXPATH_(L1:L2)//'/sex '//
     +             REFFILE(LLL1:LLL2)//',xyzxyzxyz.fits '//
     +             '-MAG_ZEROPOINT '//CDUMMY(LL1:LL2))
                ELSE
                  WRITE(*,101) '% '//SEXPATH_(L1:L2)//'/sex \\'
                  WRITE(*,100) '  xyzxyzxyz.fits -MAG_ZEROPOINT '
                  WRITE(*,101) CDUMMY(LL1:LL2)
                  ISYSTEM=SYSTEMFUNCTION(SEXPATH_(L1:L2)//'/sex '//
     +             'xyzxyzxyz.fits -MAG_ZEROPOINT '//CDUMMY(LL1:LL2))
                END IF
              ELSE
                IF(CREF.EQ.'y')THEN
                  WRITE(*,101) '% '//SEXPATH_(L1:L2)//'/sex \\'
                  WRITE(*,100) '  '//REFFILE(LLL1:LLL2)
                  WRITE(*,101) ',xyzxyzxyz.fits'
                  ISYSTEM=SYSTEMFUNCTION(SEXPATH_(L1:L2)//'/sex '//
     +             REFFILE(LLL1:LLL2)//',xyzxyzxyz.fits')
                ELSE
                  WRITE(*,100) '% '//SEXPATH_(L1:L2)//'/sex '
                  WRITE(*,101) 'xyzxyzxyz.fits'
                  ISYSTEM=SYSTEMFUNCTION(SEXPATH_(L1:L2)//'/sex '//
     +             'xyzxyzxyz.fits')
                END IF
              END IF
              ISYSTEM=SYSTEMFUNCTION('rm -f xyzxyzxyz.fits')
              WRITE(*,101) '> File test.cat has been created.'
              OPEN(44,FILE='test.cat',STATUS='OLD',FORM='FORMATTED')
10            READ(44,101,END=12) CLINEA
              IF(CLINEA(1:1).EQ.'#')THEN
                WRITE(*,101) CLINEA(1:TRUELEN(CLINEA))
              ELSE
                GOTO 12
              END IF
              GOTO 10
12            CLOSE(44)
            END IF
C..............................................................................
          ELSEIF(INDEX('i0123',COPC).NE.0)THEN
            IF(LPOSTSCRIPT)THEN
              CALL PGQID(IDOLD)
              IDNEW=PGOPEN('?')
              CALL SUBLOOK(.FALSE.,NCBUFF,.TRUE.)
            END IF
            NCOMMENTS_CAT=0
            OPEN(44,FILE='test.cat',STATUS='OLD',FORM='FORMATTED')
20          READ(44,101,END=22) CLINEA
            IF(CLINEA(1:1).EQ.'#')THEN
              WRITE(*,101) CLINEA(1:TRUELEN(CLINEA))
              NCOMMENTS_CAT=NCOMMENTS_CAT+1
              IF(CLINEA(7:22).EQ.'X_IMAGE') 
     +         READ(CLINEA(2:5),*) N_X_IMAGE
              IF(CLINEA(7:22).EQ.'Y_IMAGE')
     +         READ(CLINEA(2:5),*) N_Y_IMAGE
              IF(CLINEA(7:22).EQ.'XMIN_IMAGE') 
     +         READ(CLINEA(2:5),*) N_XMIN_IMAGE
              IF(CLINEA(7:22).EQ.'XMAX_IMAGE') 
     +         READ(CLINEA(2:5),*) N_XMAX_IMAGE
              IF(CLINEA(7:22).EQ.'YMIN_IMAGE') 
     +         READ(CLINEA(2:5),*) N_YMIN_IMAGE
              IF(CLINEA(7:22).EQ.'YMAX_IMAGE') 
     +         READ(CLINEA(2:5),*) N_YMAX_IMAGE
              IF(CLINEA(7:22).EQ.'A_IMAGE') 
     +         READ(CLINEA(2:5),*) N_A_IMAGE
              IF(CLINEA(7:22).EQ.'B_IMAGE') 
     +         READ(CLINEA(2:5),*) N_B_IMAGE
              IF(CLINEA(7:22).EQ.'THETA_IMAGE') 
     +         READ(CLINEA(2:5),*) N_THETA_IMAGE
              IF(CLINEA(7:22).EQ.'CLASS_STAR') 
     +         READ(CLINEA(2:5),*) N_CLASS_STAR
            ELSE
              GOTO 22
            END IF
            GOTO 20
22          CLOSE(44)
            IF(NCOMMENTS_CAT.GT.MAX_NCOMMENTS_CAT)THEN
              WRITE(*,100) 'NCOMMENTS_CAT....: '
              WRITE(*,*) NCOMMENTS_CAT
              WRITE(*,100) 'MAX_NCOMMENTS_CAT: '
              WRITE(*,*) MAX_NCOMMENTS_CAT
              WRITE(*,101) 'ERROR: NCOMMENTS_CAT.GT.MAX_NCOMMENTS_CAT'
              WRITE(*,100) 'Press <CR> to continue...'
              READ(*,*)
              RETURN
            END IF
C
            IF(COPC.EQ.'i') N0=READI('Object number','@')
C
            OPEN(44,FILE='test.cat',STATUS='OLD',FORM='FORMATTED')
            DO I=1,NCOMMENTS_CAT
              READ(44,101) CLINEA
              CLABEL(I)=CLINEA(1:22)
            END DO
            DO I=1,NCOMMENTS_CAT-1
              READ(CLABEL(I)(2:5),*) N1_
              READ(CLABEL(I+1)(2:5),*) N2_
              NDATA_REGISTER(I)=N2_-N1_
            END DO
            NDATA_REGISTER(NCOMMENTS_CAT)=1 !it should be always 1 (CLASS_STAR)
c
24          READ(44,101,END=26) CLINEA
            IF(COPC.EQ.'i')THEN !info
              READ(CLINEA,*) N0_
              IF(N0.EQ.N0_)THEN
                ISUM=0
                DO I=1,NCOMMENTS_CAT
                  DO II=1,NDATA_REGISTER(I)
                    ISUM=ISUM+1
                    WRITE(*,100) CLABEL(I)//': '
                    FDUMMY=LEECOLUMN(CLINEA,ISUM)
                    WRITE(*,*) FDUMMY
                  END DO
                END DO
                READ(CLINEA,*) N0
                WRITE(CDUMMY,*) N0
                L1=TRUEBEG(CDUMMY)
                L2=TRUELEN(CDUMMY)
                X_IMAGE=LEECOLUMN(CLINEA,N_X_IMAGE)
                XMIN_IMAGE=LEECOLUMN(CLINEA,N_XMIN_IMAGE)
                XMAX_IMAGE=LEECOLUMN(CLINEA,N_XMAX_IMAGE)
                YMIN_IMAGE=LEECOLUMN(CLINEA,N_YMIN_IMAGE)
                YMAX_IMAGE=LEECOLUMN(CLINEA,N_YMAX_IMAGE)
                CALL PGSCI(6)
                IF((X_IMAGE.GT.XMIN).AND.(X_IMAGE.LT.XMAX).AND.
     +           (YMAX_IMAGE.GT.YMIN).AND.(YMAX_IMAGE.LT.YMAX))THEN
                  CALL PGSLW(4)
                  CALL PGPTXT(X_IMAGE,YMAX_IMAGE,0.0,0.5,CDUMMY(L1:L2))
                  CALL PGSLW(1)
                  CALL PGMOVE(XMIN_IMAGE,YMIN_IMAGE)
                  CALL PGDRAW(XMAX_IMAGE,YMIN_IMAGE)
                  CALL PGDRAW(XMAX_IMAGE,YMAX_IMAGE)
                  CALL PGDRAW(XMIN_IMAGE,YMAX_IMAGE)
                  CALL PGDRAW(XMIN_IMAGE,YMIN_IMAGE)
                END IF
              END IF
            ELSEIF(COPC.EQ.'0')THEN !number in the catalogue
              READ(CLINEA,*) N0
              WRITE(CDUMMY,*) N0
              L1=TRUEBEG(CDUMMY)
              L2=TRUELEN(CDUMMY)
              X_IMAGE=LEECOLUMN(CLINEA,N_X_IMAGE)
              YMAX_IMAGE=LEECOLUMN(CLINEA,N_YMAX_IMAGE)
              CALL PGSCI(1)
              IF((X_IMAGE.GT.XMIN).AND.(X_IMAGE.LT.XMAX).AND.
     +         (YMAX_IMAGE.GT.YMIN).AND.(YMAX_IMAGE.LT.YMAX))THEN
                CALL PGSLW(4)
                IF(LPOSTSCRIPT)THEN
                  CALL PGQCH(CHOLD)
                  CALL PGSCH(0.7)
                END IF
                CALL PGPTXT(X_IMAGE,YMAX_IMAGE,0.0,0.5,CDUMMY(L1:L2))
                IF(LPOSTSCRIPT)THEN
                  CALL PGSCH(CHOLD)
                END IF
                CALL PGSLW(1)
              END IF
            ELSEIF(COPC.EQ.'1')THEN !barycenter
              X_IMAGE=LEECOLUMN(CLINEA,N_X_IMAGE)
              Y_IMAGE=LEECOLUMN(CLINEA,N_Y_IMAGE)
              CALL PGSCI(3)
              CALL PGPOINT(1,X_IMAGE,Y_IMAGE,2)
            ELSEIF(COPC.EQ.'2')THEN !rectangle which encloses the object
              XMIN_IMAGE=LEECOLUMN(CLINEA,N_XMIN_IMAGE)
              XMAX_IMAGE=LEECOLUMN(CLINEA,N_XMAX_IMAGE)
              YMIN_IMAGE=LEECOLUMN(CLINEA,N_YMIN_IMAGE)
              YMAX_IMAGE=LEECOLUMN(CLINEA,N_YMAX_IMAGE)
              CLASS_STAR=LEECOLUMN(CLINEA,N_CLASS_STAR)
              IF(CLASS_STAR.GT.0.8) CALL PGSLS(4)
              CALL PGSCI(4)
              CALL PGMOVE(XMIN_IMAGE,YMIN_IMAGE)
              CALL PGDRAW(XMAX_IMAGE,YMIN_IMAGE)
              CALL PGDRAW(XMAX_IMAGE,YMAX_IMAGE)
              CALL PGDRAW(XMIN_IMAGE,YMAX_IMAGE)
              CALL PGDRAW(XMIN_IMAGE,YMIN_IMAGE)
              IF(CLASS_STAR.GT.0.8) CALL PGSLS(1)
            ELSEIF(COPC.EQ.'3')THEN !isophotal limit
              X_IMAGE=LEECOLUMN(CLINEA,N_X_IMAGE)
              Y_IMAGE=LEECOLUMN(CLINEA,N_Y_IMAGE)
              A_IMAGE=LEECOLUMN(CLINEA,N_A_IMAGE)*3.0
              B_IMAGE=LEECOLUMN(CLINEA,N_B_IMAGE)*3.0
              THETA_IMAGE=-LEECOLUMN(CLINEA,N_THETA_IMAGE)*3.141593/180.
              DO I=0,360
                THETA0=REAL(I)*3.141593/180.0
                R0=A_IMAGE*B_IMAGE/SQRT(A_IMAGE*A_IMAGE-
     +           (A_IMAGE*A_IMAGE-B_IMAGE*B_IMAGE)*
     +           COS(THETA0)*COS(THETA0))
                X0=R0*COS(THETA0)
                Y0=R0*SIN(THETA0)
                XELL(I+1)=X_IMAGE
     +           +X0*COS(THETA_IMAGE)+Y0*SIN(THETA_IMAGE)
                YELL(I+1)=Y_IMAGE
     +           -X0*SIN(THETA_IMAGE)+Y0*COS(THETA_IMAGE)
              END DO
              CLASS_STAR=LEECOLUMN(CLINEA,N_CLASS_STAR)
              IF(CLASS_STAR.GT.0.8) CALL PGSLS(4)
              CALL PGSCI(5)
              CALL PGLINE(361,XELL,YELL)
              IF(CLASS_STAR.GT.0.8) CALL PGSLS(1)
            END IF
            GOTO 24
26          CLOSE(44)
            IF(LPOSTSCRIPT)THEN
              CALL PGCLOS(IDNEW)
              CALL PGSLCT(IDOLD)
            END IF
            CALL PGSCI(1)
C..............................................................................
          END IF
        END DO
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
