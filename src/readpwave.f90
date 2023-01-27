! Lee el fichero con el polinomio de calibracion en longitud de onda
        SUBROUTINE READPWAVE(INFILE,NDEGWV,CWV)
        IMPLICIT NONE
        CHARACTER*(*) INFILE
        INTEGER NDEGWV
        REAL CWV(20)
!
        INTEGER K
!------------------------------------------------------------------------------
        OPEN(40,FILE=INFILE,STATUS='OLD',FORM='FORMATTED')
        NDEGWV=0
10      READ(40,*,END=20) K,CWV(NDEGWV+1)
        NDEGWV=NDEGWV+1
        GOTO 10
20      CLOSE(40)
        NDEGWV=NDEGWV-1
!------------------------------------------------------------------------------
        END
