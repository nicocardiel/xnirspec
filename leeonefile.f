        SUBROUTINE LEEONEFILE(INFILE_,LOK)
        IMPLICIT NONE
        CHARACTER*255 INFILE_
        LOGICAL LOK
C functions
        INTEGER TRUEBEG
        INTEGER TRUELEN
        INTEGER SYSTEMFUNCTION
C variables
        INTEGER L1,L2
        INTEGER NLINES
        INTEGER ISYSTEM
        CHARACTER*255 CLINEA
        LOGICAL LOGFILE
C------------------------------------------------------------------------------
        L1=TRUEBEG(INFILE_)
        L2=TRUELEN(INFILE_)
        ISYSTEM=SYSTEMFUNCTION('ls *'//INFILE_(L1:L2)//
     +   '* > .tmp_input_file_xnirspec')
        IF(ISYSTEM.NE.0)THEN
          LOK=.FALSE.
          RETURN
        END IF
C------------------------------------------------------------------------------
        NLINES=0
10      OPEN(10,FILE='.tmp_input_file_xnirspec',STATUS='OLD',
     +   FORM='FORMATTED',ERR=99)
        READ(10,101,END=99) CLINEA
        NLINES=NLINES+1
        GOTO 10
99      CLOSE(10)
        IF(NLINES.EQ.1)THEN
          LOK=.TRUE.
          INFILE_=CLINEA
        ELSE
          LOK=.FALSE.
        END IF
        ISYSTEM=SYSTEMFUNCTION('rm -f .tmp_input_file_xnirspec')
C------------------------------------------------------------------------------
101     FORMAT(A)
        END
