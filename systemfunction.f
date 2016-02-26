        INTEGER FUNCTION SYSTEMFUNCTION(COMANDO)
        IMPLICIT NONE
        CHARACTER*(*) COMANDO
C------------------------------------------------------------------------------
        CALL SYSTEM(COMANDO,SYSTEMFUNCTION)
        END
