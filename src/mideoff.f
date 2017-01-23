C Version 23-Oct-2000
        SUBROUTINE MIDEOFF(NCBUFF)
        IMPLICIT NONE
        INTEGER NCBUFF
C
        INTEGER NBOXMAX
        PARAMETER (NBOXMAX=9)
C
        INCLUDE 'dimensions.inc'
C
        INTEGER READI,READILIM
        INTEGER TRUEBEG,TRUELEN
        CHARACTER*255 READC
C
        INTEGER L,IDUM
        INTEGER L1,L2
        INTEGER NEWBUFF1,NEWBUFF2,NEWBUFF3
        INTEGER NFRAMES(NMAXBUFF)
        INTEGER NAXISFRAME(2,9,NMAXBUFF)
        INTEGER NSIZEFB9(NMAXBUFF)
        INTEGER NF,NFRAMES_
        INTEGER NEXTRA_XMINB9,NEXTRA_XMAXB9
        INTEGER NEXTRA_YMINB9,NEXTRA_YMAXB9
        INTEGER J0(NBOXMAX),I0(NBOXMAX)
        INTEGER NXMAXB9_,NYMAXB9_
        INTEGER NPIXREGION
        INTEGER IDOLD,NPOST
        REAL FOFFSETX(NBOXMAX),FOFFSETY(NBOXMAX)
        REAL FOFFSETX_,FOFFSETY_
        REAL EFOFFSETX_,EFOFFSETY_
        REAL FJ0(NBOXMAX),FI0(NBOXMAX)
        CHARACTER*1 CPOST,CCONT
        CHARACTER*50 CBASEPOST
        CHARACTER*50 CDUMMY
        CHARACTER*80 COFFSETFILE
        LOGICAL LOGFILE
C
        COMMON/BLKNFRAMES/NFRAMES
        COMMON/BLKB9POST1/CPOST,CBASEPOST
        COMMON/BLKB9POST2/NPOST,IDOLD
        COMMON/BLKNAXISFRAME/NAXISFRAME
        COMMON/BLKNSIZEFB9/NSIZEFB9
C------------------------------------------------------------------------------
        CBASEPOST='pgplot'
C pedimos offsets
        LOGFILE=.FALSE.
        DO WHILE(.NOT.LOGFILE)
          COFFSETFILE(1:80)=
     +     READC('Name of the file with offsets (none=exit)',
     +     'offsets.dat','@')
          L1=TRUEBEG(COFFSETFILE)
          L2=TRUELEN(COFFSETFILE)
          IF(COFFSETFILE(L1:L2).EQ.'none') RETURN !permitimos escapar
          INQUIRE(FILE=COFFSETFILE(L1:L2),EXIST=LOGFILE)
          IF(.NOT.LOGFILE)THEN
            WRITE(*,101) 'ERROR: this file does not exist. Try again.'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
          END IF
        END DO
        WRITE(*,101) 'Reading offsets...'
        OPEN(20,FILE=COFFSETFILE(L1:L2),STATUS='OLD',FORM='FORMATTED')
        NF=0
5        READ(20,*,END=6) IDUM,FOFFSETX_,EFOFFSETX_,FOFFSETY_,EFOFFSETY_
        NF=NF+1
        IF(NF.LE.NBOXMAX)THEN
          FOFFSETX(NF)=FOFFSETX_
          FOFFSETY(NF)=FOFFSETY_
          WRITE(*,'(I1,2(2X,F8.3))') NF,FOFFSETX(NF),FOFFSETY(NF)
        ELSE
          CLOSE(20)
          WRITE(*,101) 'ERROR: file with offsets contains too many'//
     +     'lines (>9)'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
        GOTO 5
6       CLOSE(20)
        NFRAMES_=NF
        WRITE(*,101) '...OK! File read and closed.'
C proteccion
        IF(NFRAMES_.NE.NFRAMES(NCBUFF))THEN
          WRITE(*,100) 'NFRAMES(NCBUFF),NFRAMES_: '
          WRITE(*,*) NFRAMES(NCBUFF),NFRAMES_
          WRITE(*,100) 'ERROR: number of frames does not correspond'
          WRITE(*,101) ' with expected value.'
          CCONT(1:1)=READC('Do you want to continue anyway with '//
     +     '9x256x256 mosaic (y/n)','y','yn')
          IF(CCONT.EQ.'n') RETURN
          DO NF=1,NFRAMES_
            NAXISFRAME(1,NF,NCBUFF)=256
            NAXISFRAME(2,NF,NCBUFF)=256
          END DO
          NSIZEFB9(NCBUFF)=256
        END IF
C estimamos partes enteras y fraccionarias de los offsets, y determinamos la
C dimension maxima que vamos a necesitar para la imagen combinada
        NEXTRA_XMINB9=0
        NEXTRA_XMAXB9=0
        NEXTRA_YMINB9=0
        NEXTRA_YMAXB9=0
        DO NF=1,NFRAMES_
          J0(NF)=INT(FOFFSETX(NF))
          FJ0(NF)=FOFFSETX(NF)-REAL(J0(NF))
          IF(FOFFSETX(NF).GT.0.0)THEN
            IF(FJ0(NF).EQ.0.0)THEN
              NEXTRA_XMINB9=MIN0(NEXTRA_XMINB9,-J0(NF))
            ELSE
              NEXTRA_XMINB9=MIN0(NEXTRA_XMINB9,-J0(NF)-1)
            END IF
          ELSEIF(FOFFSETX(NF).LE.0.0)THEN
            IF(FJ0(NF).EQ.0.0)THEN
              NEXTRA_XMAXB9=MAX0(NEXTRA_XMAXB9,-J0(NF))
            ELSE
              NEXTRA_XMAXB9=MAX0(NEXTRA_XMAXB9,-J0(NF)+1)
            END IF
          END IF
          I0(NF)=INT(FOFFSETY(NF))
          FI0(NF)=FOFFSETY(NF)-REAL(I0(NF))
          IF(FOFFSETY(NF).GT.0.0)THEN
            IF(FI0(NF).EQ.0.0)THEN
              NEXTRA_YMINB9=MIN0(NEXTRA_YMINB9,-I0(NF))
            ELSE
              NEXTRA_YMINB9=MIN0(NEXTRA_YMINB9,-I0(NF)-1)
            END IF
          ELSEIF(FOFFSETY(NF).LE.0.0)THEN
            IF(FI0(NF).EQ.0.0)THEN
              NEXTRA_YMAXB9=MAX0(NEXTRA_YMAXB9,-I0(NF))
            ELSE
              NEXTRA_YMAXB9=MAX0(NEXTRA_YMAXB9,-I0(NF)+1)
            END IF
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,100) '>>> Initial extra limits of combined frame: '
        WRITE(*,*) NEXTRA_XMINB9,NEXTRA_XMAXB9,
     +   NEXTRA_YMINB9,NEXTRA_YMAXB9
        NXMAXB9_=NSIZEFB9(NCBUFF)-NEXTRA_XMINB9+NEXTRA_XMAXB9
        NYMAXB9_=NSIZEFB9(NCBUFF)-NEXTRA_YMINB9+NEXTRA_YMAXB9
C mostramos taman~o
        WRITE(*,100) '>>> New dimensions of combined frame: '
        WRITE(CDUMMY,'(I10,A1,I10)') NXMAXB9_,'x',NYMAXB9_
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101) CDUMMY(1:L)
        NPIXREGION=NXMAXB9_*NYMAXB9_
C proteccion
        IF(NXMAXB9_.GT.NXMAXB9)THEN
          WRITE(*,101) 'ERROR: X dimension > NXMAXB9'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
        IF(NYMAXB9_.GT.NYMAXB9)THEN
          WRITE(*,101) 'ERROR: Y dimension > NYMAXB9'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
C------------------------------------------------------------------------------
C definimos los buffers de trabajo
        NEWBUFF1=NCBUFF+1
        IF(NEWBUFF1.GT.NMAXBUFF/2) NEWBUFF1=1
        WRITE(CDUMMY,*) NEWBUFF1
        NEWBUFF1=READILIM(
     +   'First  auxiliary buffer # to store manipulated data..',
     +   CDUMMY,1,NMAXBUFF/2)
c
        NEWBUFF2=NEWBUFF1+1
        IF(NEWBUFF2.GT.NMAXBUFF/2) NEWBUFF2=1
        WRITE(CDUMMY,*) NEWBUFF2
        NEWBUFF2=READILIM(
     +   'Second auxiliary buffer # to store manipulated data..',
     +   CDUMMY,1,NMAXBUFF/2)
c
        NEWBUFF3=NEWBUFF2+1
        IF(NEWBUFF3.GT.NMAXBUFF/2) NEWBUFF3=1
        WRITE(CDUMMY,*) NEWBUFF3
        NEWBUFF3=READILIM(
     +   'Third auxiliary buffer # to store manipulated data...',
     +   CDUMMY,1,NMAXBUFF/2)
C------------------------------------------------------------------------------
C indicamos si queremos generar imagenes en postscript
        WRITE(*,*)
        CPOST(1:1)=READC('Create PostScript files (y/n)','n','yn')
ccc        CPOST='n'
        IF(CPOST.EQ.'y')THEN
          CBASEPOST(1:50)=READC('Base for postscript file names',
     +     CBASEPOST,'@')
          L1=TRUEBEG(CBASEPOST)
          L2=TRUELEN(CBASEPOST)
          WRITE(CDUMMY,*) NPOST+1
          NPOST=READI('Number of first'//
     +     CBASEPOST(L1:L2)//'????.ps file',CDUMMY)
          NPOST=NPOST-1
          CALL PGQID(IDOLD)
        END IF
C------------------------------------------------------------------------------
C medimos los offsets
        CALL MIDEOFFSET(NFRAMES_,NCBUFF,NEWBUFF1,NEWBUFF2,
     +   NEWBUFF3,FOFFSETX,FOFFSETY)
        CPOST='n'
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
