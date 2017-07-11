        SUBROUTINE PLOTDS9REG(NCBUFF)
        IMPLICIT NONE
        INTEGER NCBUFF
C
        INCLUDE 'dimensions.inc'
C
        INTEGER TRUELEN
C
        INTEGER L,L1,L2
        REAL X1,Y1,X2,Y2
        REAL XMIN,XMAX,YMIN,YMAX
        REAL OFFY
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
        OPEN(50,FILE=DS9REGFILE(NCBUFF),STATUS='OLD',FORM='FORMATTED')
        CALL PGBBUF
        CALL PGQWIN(XMIN,XMAX,YMIN,YMAX)
        OFFY=(YMAX-YMIN)/60.
        CALL PGSCI(NCBUFF)
10      READ(50,101,END=20) CLINEA
        IF(TRUELEN(CLINEA).GT.0)THEN
          IF(CLINEA(1:1).NE.'#')THEN
            L=INDEX(CLINEA,'line')
            IF(L.NE.0)THEN
              READ(CLINEA(L+4:),*) X1,Y1,X2,Y2
              CALL PGMOVE(X1,Y1)
              CALL PGDRAW(X2,Y2)
            ELSE
              L=INDEX(CLINEA,'text')
              IF(L.NE.0)THEN
                READ(CLINEA(L+4:),*) X1,Y1
                IF((XMIN.LT.X1).AND.(X1.LE.XMAX).AND.
     +             (YMIN.LT.Y1).AND.(Y1.LE.YMAX))THEN
                  L1=INDEX(CLINEA,'{')
                  L2=INDEX(CLINEA,'}')
                  CTEXT=CLINEA(L1+1:L2-1)
                  L=TRUELEN(CTEXT)
                  CALL PGPTEXT(X1,Y1-OFFY,0.,0.5,CTEXT(1:L))
                END IF
              END IF
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
