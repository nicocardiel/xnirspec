        INTEGER FUNCTION SYSTEMFUNCTION(COMANDO)
        IMPLICIT NONE
        CHARACTER*(*) COMANDO
!------------------------------------------------------------------------------
        CALL SYSTEM(COMANDO,SYSTEMFUNCTION)
        END
