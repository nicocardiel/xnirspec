        SUBROUTINE PLOTDS9REG(NCBUFF)
        IMPLICIT NONE
        INTEGER NCBUFF
C
        INCLUDE 'dimensions.inc'
C
        INTEGER TRUEBEG
        INTEGER TRUELEN
C
        INTEGER L,L1,L2
        INTEGER COLORDEFAULT,COLORSET
        REAL X1,Y1,X2,Y2
        REAL XMIN,XMAX,YMIN,YMAX
        REAL OFFY_TEXT
        CHARACTER*255 DS9REGFILE(NMAXBUFF)
        CHARACTER*255 CLINEA,CTEXT
        LOGICAL LOGFILE
C
        COMMON/BLKDS9REGFILE/DS9REGFILE
C------------------------------------------------------------------------------
        INQUIRE(FILE=DS9REGFILE(NCBUFF),EXIST=LOGFILE)
        IF(.NOT.LOGFILE)THEN
          WRITE(*,101) 'ERROR: the following file does not exist:'
          L=TRUELEN(DS9REGFILE(NCBUFF))
          WRITE(*,101) DS9REGFILE(NCBUFF)(1:L)
          RETURN
        END IF
C
        COLORDEFAULT=NCBUFF
C
        OPEN(50,FILE=DS9REGFILE(NCBUFF),STATUS='OLD',FORM='FORMATTED')
        READ(50,*,END=20) !skip first line
        READ(50,101,END=20) CLINEA
        L=INDEX(CLINEA,'global')
        IF(L.EQ.0)THEN
          WRITE(*,101) 'ERROR: unexpected format for ds9 region file'
        ELSE
          CALL DS9COLOR(CLINEA,COLORDEFAULT,COLORSET)
          COLORDEFAULT=COLORSET
        END IF
        CALL PGBBUF
        CALL PGQWIN(XMIN,XMAX,YMIN,YMAX)
        OFFY_TEXT=(YMAX-YMIN)/60.
10      READ(50,101,END=20) CLINEA
        IF(TRUELEN(CLINEA).GT.0)THEN !ignore blank lines
          IF(CLINEA(1:1).EQ.'#') GOTO 10
          !line
          L=INDEX(CLINEA,'line')
          IF(L.NE.0)THEN
            READ(CLINEA(L+4:),*) X1,Y1,X2,Y2
            CALL DS9COLOR(CLINEA,COLORDEFAULT,COLORSET)
            CALL PGMOVE(X1,Y1)
            CALL PGDRAW(X2,Y2)
            GOTO 10
          END IF
          !text
          L=INDEX(CLINEA,'text')
          IF(L.NE.0)THEN
            READ(CLINEA(L+4:),*) X1,Y1
            IF((XMIN.LT.X1).AND.(X1.LE.XMAX).AND.
     +         (YMIN.LT.Y1).AND.(Y1.LE.YMAX))THEN
              L1=INDEX(CLINEA,'{')
              L2=INDEX(CLINEA,'}')
              CTEXT=CLINEA(L1+1:L2-1)
              L1=TRUEBEG(CTEXT)
              L2=TRUELEN(CTEXT)
              CALL DS9COLOR(CLINEA,COLORDEFAULT,COLORSET)
              CALL PGPTEXT(X1,Y1-OFFY_TEXT,0.,0.5,CTEXT(L1:L2))
              GOTO 10
            END IF
          END IF
        END IF
        GOTO 10
20      CLOSE(50)
        CALL PGSCI(1)
        CALL PGEBUF
C
101     FORMAT(A)
        END
C
C******************************************************************************
C Auxiliary subroutine to read color
C
        SUBROUTINE DS9COLOR(CLINEA,COLORDEFAULT,COLORSET)
        IMPLICIT NONE
        CHARACTER*(*) CLINEA
        INTEGER COLORDEFAULT
        INTEGER COLORSET
C
        INTEGER TRUEBEG
        INTEGER TRUELEN
C
        INTEGER L,L1,L2
        INTEGER IR,IG,IB
        CHARACTER*255 CLINEA_
C------------------------------------------------------------------------------
        CLINEA_=CLINEA                                      !duplicate variable
        COLORSET=COLORDEFAULT
C search for 'color'
        L=INDEX(CLINEA_,'color')
        IF(L.EQ.0)THEN
          CALL PGSCI(COLORSET)
          RETURN
        END IF
        L2=TRUELEN(CLINEA_)
        IF(L2.GT.L+4)THEN
          CLINEA_=CLINEA_(L+5:)
        ELSE
          CALL PGSCI(COLORSET)
          RETURN
        END IF
C search for '='
        L=INDEX(CLINEA_,'=')
        IF(L.EQ.0)THEN
          CALL PGSCI(COLORSET)
          RETURN
        END IF
        L2=TRUELEN(CLINEA_)
        IF(L2.GT.L)THEN
          CLINEA_=CLINEA_(L+1:)
        ELSE
          CALL PGSCI(COLORSET)
          RETURN
        END IF
C search for color description
        L1=TRUEBEG(CLINEA_)
        CLINEA_=CLINEA_(L1:)
C
        IF(CLINEA_(1:1).EQ.'#')THEN
          READ(CLINEA_(2:3),'(Z2)') IR
          READ(CLINEA_(4:5),'(Z2)') IG
          READ(CLINEA_(6:7),'(Z2)') IB
          CALL PGSCR(16,REAL(IR)/255,REAL(IG)/255,REAL(IB)/255)
          COLORSET=16
        ELSEIF(CLINEA_(1:5).EQ.'white')THEN
          COLORSET=1
        ELSEIF(CLINEA_(1:3).EQ.'red')THEN
          COLORSET=2
        ELSEIF(CLINEA_(1:5).EQ.'green')THEN
          COLORSET=3
        ELSEIF(CLINEA_(1:4).EQ.'blue')THEN
          COLORSET=4
        ELSEIF(CLINEA_(1:4).EQ.'cyan')THEN
          COLORSET=5
        ELSEIF(CLINEA_(1:7).EQ.'magenta')THEN
          COLORSET=6
        ELSEIF(CLINEA_(1:6).EQ.'yellow')THEN
          COLORSET=7
        END IF
        CALL PGSCI(COLORSET)
C
        END
