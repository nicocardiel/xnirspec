        SUBROUTINE BOUNDARY(LBOUNDARY)
        IMPLICIT NONE
        LOGICAL LBOUNDARY
C
        INCLUDE 'largest.inc'
C
        INTEGER TRUEBEG,TRUELEN
        INTEGER READILIM
        INTEGER SYSTEMFUNCTION
        REAL ARCLENGTH
        REAL FPOLY
        CHARACTER*255 READC
C
        INTEGER L1,L2,LL1,LL2
        INTEGER I,L,LL,LMIN
        INTEGER ISYSTEM
        INTEGER NDEGBL(NXYMAX),NDEGBA(NXYMAX)  !polynomial degrees for boundary
        INTEGER NDEGBL00,NDEGBA00                 !maximum of NDEGBL and NDEGBA
        INTEGER NLINBL,NLINBA                  !number of arc lines in boundary
        INTEGER NREF
        REAL COEFFBL(20,NXYMAX),COEFFBA(20,NXYMAX)      !pol. coef. in boundary
        REAL XMINBL(NXYMAX),XMAXBL(NXYMAX)
        REAL YMINBL(NXYMAX),YMAXBL(NXYMAX)
        REAL XMINBA(NXYMAX),XMAXBA(NXYMAX)
        REAL YMINBA(NXYMAX),YMAXBA(NXYMAX)
        REAL XF(NXYMAX),YF(NXYMAX) !puntos a dibujar
        REAL X1,Y1,X2,Y2
        REAL XC,YC
        REAL DMIN,DIST
        REAL X0,Y0
        REAL XP(NXYMAX),YP(NXYMAX)
        CHARACTER*1 COPC,CH,CSURE
        CHARACTER*50 CDUMX,CDUMY
        CHARACTER*255 FILENAME
        LOGICAL LOGFILE,LOGFILERR
        LOGICAL LBOUNDL,LBOUNDA
C
        COMMON/BLKBOUND1/NLINBL,NDEGBL,NDEGBL00
        COMMON/BLKBOUND1B/NLINBA,NDEGBA,NDEGBA00
        COMMON/BLKBOUND2/COEFFBL
        COMMON/BLKBOUND2B/COEFFBA
        COMMON/BLKBOUND3/LBOUNDL,LBOUNDA
        COMMON/BLKBOUND4/XMINBL,YMINBL,XMAXBL,YMAXBL
        COMMON/BLKBOUND5/XMINBA,YMINBA,XMAXBA,YMAXBA
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,101) '(l) load boundary from file'
        IF((NLINBA.EQ.0).AND.(NLINBL.GE.2))THEN
          WRITE(*,101) '(a) define vertical pattern of BA lines'
        END IF
        IF(LBOUNDARY)THEN
          WRITE(*,101) '(s) save boundary into file'
          WRITE(*,101) '(r) delete BL line with mouse'
          WRITE(*,101) '(d) delete BA line with mouse'
          WRITE(*,101) '(1) test shrinking/stretching in BL'
          WRITE(*,101) '(2) test shrinking/stretching in BA'
        END IF
        WRITE(*,101) '(x) exit'
        IF(LBOUNDARY)THEN
          COPC(1:1)=READC('Option','x','lsrd12x')
        ELSE
          IF((NLINBA.EQ.0).AND.(NLINBL.GE.2))THEN
            COPC(1:1)=READC('Option','x','lax')
          ELSE
            COPC(1:1)=READC('Option','x','lx')
          END IF
        END IF
C------------------------------------------------------------------------------
        IF(COPC.EQ.'x')THEN
          RETURN
C------------------------------------------------------------------------------
C define vertical pattern of BA lines
        ELSEIF(COPC.EQ.'a')THEN
          NLINBA=READILIM('No. of BA lines to be defined','@',2,100)
          NDEGBA00=READILIM('Polynomial degree for BA lines','@',0,19)
          DO L=1,NLINBA
            NDEGBA(L)=NDEGBA00
            COEFFBA(1,L)=XMINBL(1)+REAL(L)*(XMAXBL(1)-XMINBL(1))/
     +       REAL(NLINBA+1)
            IF(NDEGBA00.GE.1)THEN
              DO I=2,NDEGBA00+1
                COEFFBA(I,L)=0.00
              END DO
            END IF
            LBOUNDARY=.TRUE.
            LBOUNDL=.TRUE.
            LBOUNDA=.TRUE.
          END DO
          CALL SORTLINES
C------------------------------------------------------------------------------
C load boundary from file
        ELSEIF(COPC.EQ.'l')THEN
          LOGFILE=.FALSE.
          LOGFILERR=.FALSE.
          DO WHILE(.NOT.LOGFILE)
            FILENAME=READC('Boundary file name (none=EXIT)',
     +       '*.boun','@')
            IF((INDEX(FILENAME,'*').NE.0).OR.
     +       (INDEX(FILENAME,'?').NE.0))THEN
              L1=TRUEBEG(FILENAME)
              L2=TRUELEN(FILENAME)
              ISYSTEM=SYSTEMFUNCTION('ls '//FILENAME(L1:L2))
            ELSEIF(FILENAME.EQ.'none')THEN
              LOGFILE=.TRUE.
              LOGFILERR=.TRUE.
            ELSE
              INQUIRE(FILE=FILENAME,EXIST=LOGFILE)
              IF(.NOT.LOGFILE)THEN
                WRITE(*,101) '***ERROR***'
                WRITE(*,100) '=> This file does not exist.'
                WRITE(*,100) ' Try again. (press <CR>...)'
                READ(*,*)
              END IF
            END IF
          END DO
C
          IF(.NOT.LOGFILERR)THEN                   !the file name is not 'none'
            WRITE(*,100) 'Reading file...'
            OPEN(10,FILE=FILENAME,STATUS='OLD',FORM='FORMATTED')
c
            READ(10,*) !skip line
            READ(10,*) NLINBL
            DO L=1,NLINBL
              READ(10,*) NDEGBL(L)
              DO I=1,NDEGBL(L)+1
                READ(10,*) COEFFBL(I,L)
              END DO
              READ(10,*) XMINBL(L),YMINBL(L),XMAXBL(L),YMAXBL(L)
            END DO
c
            READ(10,*) !skip line
            READ(10,*) NLINBA
            DO L=1,NLINBA
              READ(10,*) NDEGBA(L)
              DO I=1,NDEGBA(L)+1
                READ(10,*) COEFFBA(I,L)
              END DO
              READ(10,*) XMINBA(L),YMINBA(L),XMAXBA(L),YMAXBA(L)
            END DO
c
            NDEGBL00=0
            DO L=1,NLINBL
              IF(NDEGBL(L).GT.NDEGBL00) NDEGBL00=NDEGBL(L)
            END DO
            NDEGBA00=0
            DO L=1,NLINBA
              IF(NDEGBA(L).GT.NDEGBA00) NDEGBA00=NDEGBA(L)
            END DO
c
            CLOSE(10)
            WRITE(*,101) ' ...OK!'
            LBOUNDARY=.TRUE.
            LBOUNDL=.TRUE.
            LBOUNDA=.TRUE.
          END IF
C------------------------------------------------------------------------------
C save boundary into file
        ELSEIF(COPC.EQ.'s')THEN
          LOGFILE=.TRUE.
          LOGFILERR=.FALSE.
          DO WHILE(LOGFILE)
            FILENAME=READC('Boundary file name (none=EXIT)',
     +       '*.boun','@')
            IF((INDEX(FILENAME,'*').NE.0).OR.
     +       (INDEX(FILENAME,'?').NE.0))THEN
              L1=TRUEBEG(FILENAME)
              L2=TRUELEN(FILENAME)
              ISYSTEM=SYSTEMFUNCTION('ls '//FILENAME(L1:L2))
            ELSEIF(FILENAME.EQ.'none')THEN
              LOGFILE=.FALSE.
              LOGFILERR=.TRUE.
            ELSE
              INQUIRE(FILE=FILENAME,EXIST=LOGFILE)
              IF(LOGFILE)THEN
                WRITE(*,101) '***ERROR***'
                WRITE(*,100) '=> This file does already exist.'
                WRITE(*,100) ' Try again. (press <CR>...)'
                READ(*,*)
              END IF
            END IF
          END DO
C
          IF(.NOT.LOGFILERR)THEN                   !the file name is not 'none'
            WRITE(*,100) 'Writting file...'
            OPEN(10,FILE=FILENAME,STATUS='NEW',FORM='FORMATTED')
c
            WRITE(10,101) '# Boundary: BL lines'
            WRITE(10,*) NLINBL
            DO L=1,NLINBL
              WRITE(10,*) NDEGBL(L)
              DO I=1,NDEGBL(L)+1
                WRITE(10,*) COEFFBL(I,L)
              END DO
              WRITE(10,*) XMINBL(L),YMINBL(L),XMAXBL(L),YMAXBL(L)
            END DO
c
            WRITE(10,101) '# Boundary: BA lines'
            WRITE(10,*) NLINBA
            DO L=1,NLINBA
              WRITE(10,*) NDEGBA(L)
              DO I=1,NDEGBA(L)+1
                WRITE(10,*) COEFFBA(I,L)
              END DO
              WRITE(10,*) XMINBA(L),YMINBA(L),XMAXBA(L),YMAXBA(L)
            END DO
c
            CLOSE(10)
            WRITE(*,101) ' ...OK!'
          END IF
C------------------------------------------------------------------------------
C delete BL line with mouse
        ELSEIF(COPC.EQ.'r')THEN
          IF(NLINBL.EQ.2)THEN
            WRITE(*,100) 'WARNING: only 2 BL lines defined.'
            WRITE(*,101) ' You cannot remove any BL line.'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
            RETURN
          END IF
          WRITE(*,101) 'Press mouse button in BL line to be removed...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          IF(CH.EQ.'A')THEN
            DMIN=10.*REAL(NXYMAX)
            LMIN=0
            DO L=1,NLINBL
              CALL DISTPP(XC,YC,NDEGBL(L),COEFFBL(1,L),XC,DIST,X0,Y0)
              IF(DIST.LT.DMIN)THEN
                DMIN=DIST
                LMIN=L
              END IF
            END DO
            IF(LMIN.EQ.0)THEN
              STOP 'FATAL ERROR: in subroutine BOUNDARY' !no debe pasar nunca
            ELSE
              DO I=1,NXYMAX
                XP(I)=XMINBL(LMIN)+REAL(I-1)/REAL(NXYMAX-1)*
     +           (XMAXBL(LMIN)-XMINBL(LMIN))
                YP(I)=FPOLY(NDEGBL(LMIN),COEFFBL(1,LMIN),XP(I))
              END DO
              CALL PGSCI(7)
              CALL PGLINE(NXYMAX,XP,YP)
              CALL PGSCI(1)
              CSURE(1:1)=READC(
     +         'Do you really want to eliminate this BL line','n','yn')
              IF(CSURE.EQ.'y')THEN
                IF(LMIN.LT.NLINBL)THEN
                  DO L=LMIN,NLINBL-1
                    NDEGBL(L)=NDEGBL(L+1)
                    DO I=1,NDEGBL(L+1)+1
                      COEFFBL(I,L)=COEFFBL(I,L+1)
                    END DO
                    XMINBL(L)=XMINBL(L+1)
                    XMAXBL(L)=XMAXBL(L+1)
                    YMINBL(L)=YMINBL(L+1)
                    YMAXBL(L)=YMAXBL(L+1)
                  END DO
                END IF
                NLINBL=NLINBL-1
              END IF
            END IF
          ELSE
            WRITE(*,101) '...removal canceled!'
          END IF
C------------------------------------------------------------------------------
C delete BA line with mouse
        ELSEIF(COPC.EQ.'d')THEN
          IF(NLINBA.EQ.2)THEN
            WRITE(*,100) 'WARNING: only 2 BA lines defined.'
            WRITE(*,101) ' You cannot remove any BA line.'
            WRITE(*,100) 'Press <CR> to continue...'
            READ(*,*)
            RETURN
          END IF
          WRITE(*,101) 'Press mouse button in BA line to be removed...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          IF(CH.EQ.'A')THEN
            DMIN=10.*REAL(NXYMAX)
            LMIN=0
            DO L=1,NLINBA
              CALL DISTPP(YC,XC,NDEGBA(L),COEFFBA(1,L),YC,DIST,Y0,X0)
              IF(DIST.LT.DMIN)THEN
                DMIN=DIST
                LMIN=L
              END IF
            END DO
            IF(LMIN.EQ.0)THEN
              STOP 'FATAL ERROR: in subroutine BOUNDARY' !no debe pasar nunca
            ELSE
              DO I=1,NXYMAX
                YP(I)=YMINBA(LMIN)+REAL(I-1)/REAL(NXYMAX-1)*
     +           (YMAXBA(LMIN)-YMINBA(LMIN))
                XP(I)=FPOLY(NDEGBA(LMIN),COEFFBA(1,LMIN),YP(I))
              END DO
              CALL PGSCI(7)
              CALL PGLINE(NXYMAX,XP,YP)
              CALL PGSCI(1)
              CSURE(1:1)=READC(
     +         'Do you really want to eliminate this BA line','n','yn')
              IF(CSURE.EQ.'y')THEN
                IF(LMIN.LT.NLINBA)THEN
                  DO L=LMIN,NLINBA-1
                    NDEGBA(L)=NDEGBA(L+1)
                    DO I=1,NDEGBA(L+1)+1
                      COEFFBA(I,L)=COEFFBA(I,L+1)
                    END DO
                    XMINBA(L)=XMINBA(L+1)
                    XMAXBA(L)=XMAXBA(L+1)
                    YMINBA(L)=YMINBA(L+1)
                    YMAXBA(L)=YMAXBA(L+1)
                  END DO
                END IF
                NLINBA=NLINBA-1
              END IF
            END IF
          ELSE
            WRITE(*,101) '...removal canceled!'
          END IF
C------------------------------------------------------------------------------
C test shrinking/stretching in BL
        ELSEIF(COPC.EQ.'1')THEN
          NREF=READILIM('Reference BL line','1',1,NLINBL)
          WRITE(CDUMX,*) NREF
          L1=TRUEBEG(CDUMX)
          L2=TRUELEN(CDUMX)
C XF va a ser la escala de referencia
C YF va a ser la escala que vamos a ajustar frente a XF
          CALL INTERSEC(NDEGBA(1),COEFFBA(1,1),
     +     NDEGBL(NREF),COEFFBL(1,NREF),
     +     0.5*(XMINBA(1)+XMAXBA(1)),X1,Y1)
          XF(1)=0.
          DO L=2,NLINBA
            CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),
     +       NDEGBL(NREF),COEFFBL(1,NREF),
     +       0.5*(XMINBA(L)+XMAXBA(L)),X2,Y2)
            XF(L)=ARCLENGTH(NDEGBL(NREF),COEFFBL(1,NREF),
     +       X1,X2,100)
          END DO
          DO LL=1,NLINBL
            IF(LL.NE.NREF)THEN !si no es la escala de referencia la ajustamos
              WRITE(CDUMY,*) LL
              LL1=TRUEBEG(CDUMY)
              LL2=TRUELEN(CDUMY)
              CALL INTERSEC(NDEGBA(1),COEFFBA(1,1),
     +         NDEGBL(LL),COEFFBL(1,LL),
     +         0.5*(XMINBA(1)+XMAXBA(1)),X1,Y1)
              YF(1)=0.
              DO L=2,NLINBA
                CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),
     +           NDEGBL(LL),COEFFBL(1,LL),
     +           0.5*(XMINBA(L)+XMAXBA(L)),X2,Y2)
                YF(L)=ARCLENGTH(NDEGBL(LL),COEFFBL(1,LL),
     +           X1,X2,100)-XF(L)
              END DO
              CALL SUBPLOT(NLINBA,1,NLINBA,XF,YF,XF,YF,
     +         .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +         'reference scale BL#'//CDUMX(L1:L2),
     +         '(measured-reference) scale BL#'//CDUMY(LL1:LL2),
     +         'test shrinking/stretching',3,17,1.5)
              WRITE(*,*)
              WRITE(*,101) '=> Reference is BL#'//CDUMX(L1:L2)
              WRITE(*,101) '=> Working with BL#'//CDUMY(LL1:LL2)
              CALL SFITPOL(NLINBA,XF,YF)
            END IF
          END DO
C------------------------------------------------------------------------------
C test shrinking/stretching in BA
        ELSEIF(COPC.EQ.'2')THEN
          NREF=READILIM('Reference BA line','1',1,NLINBA)
          WRITE(CDUMX,*) NREF
          L1=TRUEBEG(CDUMX)
          L2=TRUELEN(CDUMX)
C XF va a ser la escala de referencia
C YF va a ser la escala que vamos a ajustar frente a XF
          CALL INTERSEC(NDEGBA(NREF),COEFFBA(1,NREF),
     +     NDEGBL(1),COEFFBL(1,1),
     +     0.5*(XMINBA(NREF)+XMAXBA(NREF)),X1,Y1)
          XF(1)=0.
          DO L=2,NLINBL
            CALL INTERSEC(NDEGBA(NREF),COEFFBA(1,NREF),
     +       NDEGBL(L),COEFFBL(1,L),
     +       0.5*(XMINBA(NREF)+XMAXBA(NREF)),X2,Y2)
            XF(L)=ARCLENGTH(NDEGBA(NREF),COEFFBA(1,NREF),
     +       Y1,Y2,100)
          END DO
          DO LL=1,NLINBA
            IF(LL.NE.NREF)THEN !si no es la escala de referencia la ajustamos
              WRITE(CDUMY,*) LL
              LL1=TRUEBEG(CDUMY)
              LL2=TRUELEN(CDUMY)
              CALL INTERSEC(NDEGBA(LL),COEFFBA(1,LL),
     +         NDEGBL(1),COEFFBL(1,1),
     +         0.5*(XMINBA(LL)+XMAXBA(LL)),X1,Y1)
              YF(1)=0.
              DO L=2,NLINBL
                CALL INTERSEC(NDEGBA(LL),COEFFBA(1,LL),
     +           NDEGBL(L),COEFFBL(1,L),
     +           0.5*(XMINBA(LL)+XMAXBA(LL)),X2,Y2)
                YF(L)=ARCLENGTH(NDEGBA(LL),COEFFBA(1,LL),
     +           Y1,Y2,100)-XF(L)
              END DO
              CALL SUBPLOT(NLINBL,1,NLINBL,XF,YF,XF,YF,
     +         .TRUE.,.TRUE.,.FALSE.,.FALSE.,
     +         'reference scale BA#'//CDUMX(L1:L2),
     +         '(measured-reference) scale BA#'//CDUMY(LL1:LL2),
     +         'test shrinking/stretching',3,17,1.5)
              WRITE(*,*)
              WRITE(*,101) '=> Reference is BA#'//CDUMX(L1:L2)
              WRITE(*,101) '=> Working with BA#'//CDUMY(LL1:LL2)
              CALL SFITPOL(NLINBL,XF,YF)
            END IF
          END DO
C------------------------------------------------------------------------------
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
