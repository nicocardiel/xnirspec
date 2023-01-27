! Si NSEED=0, la subrutina retorna el siguiente numero aleatorio de la
! lista. Si NSEED es negativo, entonces se aleatoriza la lista de
! numeros aleatorios utilizando el tiempo del ordenador (la semilla no
! cambia en escalas inferiores a un segundo de tiempo). Si la semilla es
! positiva, se utiliza dicha semilla para iniciar una nueva secuencia de 
! numeros aleatorios. En los dos ultimos casos NSEED se modifica y se hace 
! igual a cero para que la siguientes llamadas continuen con los numeros de 
! la nueva lista).
        REAL FUNCTION RANDOMNUMBER(NSEED)
        IMPLICIT NONE
        INTEGER NSEED
!------------------------------------------------------------------------------
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
!
        END
