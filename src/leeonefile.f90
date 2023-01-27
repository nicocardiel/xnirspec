        SUBROUTINE LEEONEFILE(INFILE_,LOK)
        IMPLICIT NONE
        CHARACTER*255 INFILE_
        LOGICAL LOK
! functions
        INTEGER TRUEBEG
        INTEGER TRUELEN
        INTEGER SYSTEMFUNCTION
! variables
        INTEGER L1,L2
        INTEGER NLINES
        INTEGER ISYSTEM
        CHARACTER*255 CLINEA
!------------------------------------------------------------------------------
        L1=TRUEBEG(INFILE_)
        L2=TRUELEN(INFILE_)
        ISYSTEM=SYSTEMFUNCTION('ls -l *'//INFILE_(L1:L2)//'* > .tmp_input_file_xnirspec')
        IF(ISYSTEM.NE.0)THEN
          LOK=.FALSE.
          RETURN
        END IF
!------------------------------------------------------------------------------
        NLINES=0
10      OPEN(10,FILE='.tmp_input_file_xnirspec',STATUS='OLD',FORM='FORMATTED',ERR=98)
        READ(10,101,END=98) CLINEA
        NLINES=NLINES+1
        GOTO 10
98      CLOSE(10)
        IF(NLINES.NE.1)THEN
          ISYSTEM=SYSTEMFUNCTION('cat .tmp_input_file_xnirspec')
        END IF
        ISYSTEM=SYSTEMFUNCTION('rm -f .tmp_input_file_xnirspec')
!
        IF(NLINES.NE.1)THEN
          LOK=.FALSE.
          RETURN
        END IF
!
        ISYSTEM=SYSTEMFUNCTION('ls *'//INFILE_(L1:L2)//'* > .tmp_input_file_xnirspec_single')
        OPEN(20,FILE='.tmp_input_file_xnirspec_single',STATUS='OLD',FORM='FORMATTED',ERR=99)
        READ(20,101,END=99) CLINEA
        CLOSE(20)
        ISYSTEM=SYSTEMFUNCTION('rm -f .tmp_input_file_xnirspec_single')
        INFILE_=CLINEA
        LOK=.TRUE.
        RETURN
!
99      LOK=.FALSE.
        RETURN
!------------------------------------------------------------------------------
101     FORMAT(A)
        END
