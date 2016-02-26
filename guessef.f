C------------------------------------------------------------------------------
        SUBROUTINE GUESSEF(INFILE,OUTFILE)
        IMPLICIT NONE
        CHARACTER*(*) INFILE
        CHARACTER*(*) OUTFILE
C
        INTEGER TRUELEN
        INTEGER LIN,LOUT,I
C------------------------------------------------------------------------------
        LIN=TRUELEN(INFILE)
        LOUT=LEN(OUTFILE)
        IF(LIN.GT.LOUT-1) LIN=LOUT-1
        IF(LIN.GT.255) LIN=254
C
        DO I=1,LOUT
          OUTFILE(I:I)=' '
        END DO
C
        IF(INDEX(INFILE(1:LIN),'.').EQ.0)THEN
          OUTFILE(1:LIN)=INFILE(1:LIN)
          OUTFILE(LIN+1:LIN+1)='e'
          RETURN
        ENDIF
C
        DO I=LIN,1,-1
          IF(INFILE(I:I).EQ.'.')GOTO 10
          IF(INFILE(I:I).EQ.'/')THEN     !si hay barra de directorio, escapamos
            OUTFILE(1:LIN)=INFILE(1:LIN)
            OUTFILE(LIN+1:LIN+1)='e'
            RETURN
          END IF
        END DO
10      OUTFILE(1:I-1)=INFILE(1:I-1)
        OUTFILE(I:I)='e'
        OUTFILE(I+1:LIN+1)=INFILE(I:LIN)
C
        END
