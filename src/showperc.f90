!------------------------------------------------------------------------------
! Version 28-July-1998                                         File: showperc.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!omment
!
! SUBROUTINE SHOWPERC(N1,N2,ISTEP,I,NEXTINFO)
!
! Input: N1,N2,ISTEP,I,NEXTINFO
! Output: NEXTINFO
!
! Display the percentage of work performed in a loop, which has been defined
! in the range from N1 to N2, with an incremental step ISTEP, being I the 
! current value of the loop control.
!
! INTEGER N1 -> first limit of the control variable of the loop
! INTEGER N2 -> last limit of the control variable of the loop
! INTEGER ISTEP -> the value by which the control variable is incremented
! INTEGER I -> current value of the control variable
! INTEGER NEXTINFO -> integer which stores the last fraction of 10 displayed;
!                     this variable is initialized in the first call of this
!                     routine
!
!omment
!------------------------------------------------------------------------------
        SUBROUTINE SHOWPERC(N1,N2,ISTEP,I,NEXTINFO)
        IMPLICIT NONE
        INTEGER N1,N2,ISTEP,I,NEXTINFO
!
        INTEGER SYSTEMFUNCTION
!
        INTEGER ISYSTEM
        REAL FN1,FN2,FISTEP,FNEXTINFO
!------------------------------------------------------------------------------
        IF(NEXTINFO.GT.10) RETURN
!
        IF(I.EQ.N1)THEN
          WRITE(*,'(A)') 'Working...'
          NEXTINFO=0
        ELSE
          IF(I.LT.N1) RETURN
        END IF
!
        FN1=REAL(N1)
        FN2=REAL(N2)
        FISTEP=REAL(ISTEP)
        FNEXTINFO=REAL(NEXTINFO)
!
        IF((FN2-FN1)/FISTEP.GT.10.)THEN
          IF(I-N1.GE.NINT((FN2-FN1)/10.*FNEXTINFO))THEN
            IF(NEXTINFO.EQ.10)THEN
              WRITE(*,'(A8)') '100% OK!'
              NEXTINFO=NEXTINFO+1
            ELSE
              WRITE(*,'(I3,A6,$)')10*NEXTINFO,'%,... '
              ISYSTEM=SYSTEMFUNCTION('date')
              NEXTINFO=NEXTINFO+1
            END IF
          END IF
        ELSE
          IF(I.EQ.N2)THEN
            WRITE(*,'(I3,A4)') I,' OK!'
          ELSE
            WRITE(*,'(I3,A4)') I,',...'
          END IF
        END IF
!
        END
