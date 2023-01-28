        SUBROUTINE XCUT(NCBUFF,IY1,IY2,SX)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
!
        INCLUDE 'dimensions.inc'
        INTEGER NCBUFF
        INTEGER IY1,IY2
        REAL SX(NXMAX)
!
        INTEGER I1,I2
        INTEGER I,J
        INTEGER NAXIS(2,NMAXBUFF)
        REAL F
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
!
        COMMON/BLKNAXIS/NAXIS
!delete COMMON/BLKIMAGEN1/IMAGEN
!       COMMON/BLKIMAGEN2/NCBUFF
!------------------------------------------------------------------------------
        IF(IY1.LE.IY2)THEN
          I1=IY1
          I2=IY2
        ELSE
          I1=IY2
          I2=IY1
        END IF
!
        DO J=1,NAXIS(1,NCBUFF)
          SX(J)=0.
        END DO
        DO I=I1,I2
          DO J=1,NAXIS(1,NCBUFF)
            SX(J)=SX(J)+IMAGEN(J,I,NCBUFF)
          END DO
        END DO
!
        IF(I1.EQ.I2) RETURN
!
        F=1./REAL(I2-I1+1)
        DO J=1,NAXIS(1,NCBUFF)
          SX(J)=SX(J)*F
        END DO
!
        END
