        SUBROUTINE HPLOT(COMMAND)
        IMPLICIT NONE
        CHARACTER*(*) COMMAND
C
        INCLUDE 'largest.inc'
C
        INTEGER TRUEBEG,TRUELEN
C
        INTEGER L1,L2
        INTEGER NLAST,N1LAST,N2LAST
        INTEGER ICOLORLAST,NSYMBLAST
        REAL XPLAST(NXYMAX*NOVERSAMPMAX),YPLAST(NXYMAX*NOVERSAMPMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XV1,XV2,YV1,YV2
        REAL XW1,XW2,YW1,YW2
        REAL OLD_CH
        REAL XC,YC,XC1,XC2,YC1,YC2
        REAL CHLAST
        CHARACTER*1 CH
        CHARACTER*255 CLABX,CLABY,CLABG
C
        COMMON/BLKPLIMITS/XMIN,XMAX,YMIN,YMAX
        COMMON/BLKPNLAST/NLAST,N1LAST,N2LAST
        COMMON/BLKPLAST/XPLAST,YPLAST
        COMMON/BLKPCOLSYMB/ICOLORLAST,NSYMBLAST
        COMMON/BLKPCH/CHLAST
        COMMON/BLKPLABELS/CLABX,CLABY,CLABG
C------------------------------------------------------------------------------
C comenzamos buffering
        CALL PGBBUF
C almacenamos region de dibujo actual
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(XW1,XW2,YW1,YW2)
        CALL PGQCH(OLD_CH)
C definimos nueva region de dibujo para el plot
        CALL PGSCH(0.8)
        CALL PGSVP(0.05,0.37,0.46,0.75)
C definimos las nuevas coordenadas
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
C------------------------------------------------------------------------------
        IF(COMMAND.EQ.'zoom')THEN
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          XC1=XC
          YC1=YC
          CALL PGSCI(5)
          CALL RPGBAND(2,0,XC1,YC1,XC,YC,CH)
          CALL PGSCI(1)
          XC2=XC
          YC2=YC
          IF(XC2.LT.XC1)THEN
            XC=XC2
            XC2=XC1
            XC1=XC
          END IF
          IF(YC2.LT.YC1)THEN
            YC=YC2
            YC2=YC1
            YC1=YC
          END IF
          IF(XC1.LT.XMIN) XC1=XMIN
          IF(YC1.LT.YMIN) YC1=YMIN
          IF(XC2.GT.XMAX) XC2=XMAX
          IF(YC2.GT.YMAX) YC2=YMAX
          XMIN=XC1
          XMAX=XC2
          YMIN=YC1
          YMAX=YC2
          N1LAST=1
          N2LAST=NLAST
C..............................................................................
        ELSEIF(COMMAND.EQ.'whole')THEN
          CALL FINDMML(NLAST,1,NLAST,XPLAST,XMIN,XMAX)
          CALL FINDMML(NLAST,1,NLAST,YPLAST,YMIN,YMAX)
          DX=XMAX-XMIN
          XMIN=XMIN-DX/20.
          XMAX=XMAX+DX/20.
          DY=YMAX-YMIN
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
          N1LAST=1
          N2LAST=NLAST
C..............................................................................
        ELSEIF(COMMAND.EQ.'min,max')THEN
          CALL FINDMML(NLAST,N1LAST,N2LAST,YPLAST,YMIN,YMAX)
          DY=YMAX-YMIN
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
C..............................................................................
        ELSE
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> unexpected COMMAND variable in HPLOT'
          WRITE(*,100) '=> COMMAND='
          L1=TRUEBEG(COMMAND)
          L2=TRUELEN(COMMAND)
          WRITE(*,101) COMMAND(L1:L2)
          STOP
        END IF
C------------------------------------------------------------------------------
C borramos plot previo
        CALL RPGERASW(0.00,0.40,0.38,0.80,0)
C dibujamos caja y datos
        CALL PGSWIN(XMIN,XMAX,YMIN,YMAX)
        CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
        CALL PGLABEL(CLABX,CLABY,' ')
        CALL PGMTXT('T',1.0,0.5,0.5,CLABG)
        CALL PGSCI(ICOLORLAST)
        CALL PGSCH(CHLAST)
        IF(NSYMBLAST.LE.100)THEN
          CALL PGPOINT(N2LAST-N1LAST+1,XPLAST(N1LAST),YPLAST(N1LAST),
     +     NSYMBLAST)
        ELSEIF(NSYMBLAST.LE.105)THEN
          CALL PGSLS(NSYMBLAST-100)
          CALL PGLINE(N2LAST-N1LAST+1,XPLAST(N1LAST),YPLAST(N1LAST))
          CALL PGSLS(1)
        ELSEIF(NSYMBLAST.LE.205)THEN
          CALL PGSLS(NSYMBLAST-200)
          CALL PGBIN(N2LAST-N1LAST+1,XPLAST(N1LAST),YPLAST(N1LAST),
     +     .TRUE.)
          CALL PGSLS(1)
        END IF
        CALL PGSCI(1)
C------------------------------------------------------------------------------
C recuperamos region de dibujo inicial
        CALL PGSVP(XV1,XV2,YV1,YV2)
        CALL PGSWIN(XW1,XW2,YW1,YW2)
        CALL PGSCH(OLD_CH)
C terminamos buffering
        CALL PGEBUF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
