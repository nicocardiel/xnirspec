        SUBROUTINE YCUT(NCBUFF,IX1,IX2,SY)
        USE Dynamic_Array_IMAGEN
        IMPLICIT NONE
        INCLUDE 'interface_imagen.inc'
!
        INCLUDE 'dimensions.inc'
        INTEGER NCBUFF
        INTEGER IX1,IX2
        REAL SY(NYMAX)
!
        INTEGER J1,J2
        INTEGER I,J
        INTEGER NAXIS(2,NMAXBUFF)
        REAL F
!delete REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
!
        COMMON/BLKNAXIS/NAXIS
!delete COMMON/BLKIMAGEN1/IMAGEN
!       COMMON/BLKIMAGEN2/NCBUFF
!------------------------------------------------------------------------------
        IF(IX1.LE.IX2)THEN
          J1=IX1
          J2=IX2
        ELSE
          J1=IX2
          J2=IX1
        END IF
!
        DO I=1,NAXIS(2,NCBUFF)
          SY(I)=0.
        END DO
        DO I=1,NAXIS(2,NCBUFF)
          DO J=J1,J2
            SY(I)=SY(I)+IMAGEN(J,I,NCBUFF)
          END DO
        END DO
!
        IF(J1.EQ.J2) RETURN
!
        F=1./REAL(J2-J1+1)
        DO I=1,NAXIS(2,NCBUFF)
          SY(I)=SY(I)*F
        END DO
!
        END
