C Si NSEED=0, la subrutina retorna el siguiente numero aleatorio de la
C lista. Si NSEED es negativo, entonces se aleatoriza la lista de
C numeros aleatorios utilizando el tiempo del ordenador (la semilla no
C cambia en escalas inferiores a un segundo de tiempo). Si la semilla es
C positiva, se utiliza dicha semilla para iniciar una nueva secuencia de 
C numeros aleatorios. En los dos ultimos casos NSEED se modifica y se hace 
C igual a cero para que la siguientes llamadas continuen con los numeros de 
C la nueva lista).
        REAL FUNCTION RANDOMNUMBER(NSEED)
        IMPLICIT NONE
        INTEGER NSEED
C------------------------------------------------------------------------------
        IF(NSEED.EQ.0)THEN
          RANDOMNUMBER=RAND()
          RETURN
        ELSEIF(NSEED.LT.0)THEN
          CALL SRAND(TIME())
        ELSEIF(NSEED.GT.0)THEN
          CALL SRAND(NSEED)
        END IF
        NSEED=0
        RANDOMNUMBER=RAND()
C
        END
