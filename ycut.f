        SUBROUTINE YCUT(NCBUFF,IX1,IX2,SY)
        IMPLICIT NONE
        INCLUDE 'dimensions.inc'
        INTEGER NCBUFF
        INTEGER IX1,IX2
        REAL SY(NYMAX)
C
        INTEGER J1,J2
        INTEGER I,J
        INTEGER NAXIS(2,NMAXBUFF)
        REAL F
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
C
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKIMAGEN1/IMAGEN
!       COMMON/BLKIMAGEN2/NCBUFF
C------------------------------------------------------------------------------
        IF(IX1.LE.IX2)THEN
          J1=IX1
          J2=IX2
        ELSE
          J1=IX2
          J2=IX1
        END IF
C
        DO I=1,NAXIS(2,NCBUFF)
          SY(I)=0.
        END DO
        DO I=1,NAXIS(2,NCBUFF)
          DO J=J1,J2
            SY(I)=SY(I)+IMAGEN(J,I,NCBUFF)
          END DO
        END DO
C
        IF(J1.EQ.J2) RETURN
C
        F=1./REAL(J2-J1+1)
        DO I=1,NAXIS(2,NCBUFF)
          SY(I)=SY(I)*F
        END DO
C
        END
