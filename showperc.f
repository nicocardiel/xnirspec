C------------------------------------------------------------------------------
C Version 28-July-1998                                         File: showperc.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE SHOWPERC(N1,N2,ISTEP,I,NEXTINFO)
C
C Input: N1,N2,ISTEP,I,NEXTINFO
C Output: NEXTINFO
C
C Display the percentage of work performed in a loop, which has been defined
C in the range from N1 to N2, with an incremental step ISTEP, being I the 
C current value of the loop control.
C
C INTEGER N1 -> first limit of the control variable of the loop
C INTEGER N2 -> last limit of the control variable of the loop
C INTEGER ISTEP -> the value by which the control variable is incremented
C INTEGER I -> current value of the control variable
C INTEGER NEXTINFO -> integer which stores the last fraction of 10 displayed;
C                     this variable is initialized in the first call of this
C                     routine
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE SHOWPERC(N1,N2,ISTEP,I,NEXTINFO)
        IMPLICIT NONE
        INTEGER N1,N2,ISTEP,I,NEXTINFO
C
        REAL FN1,FN2,FISTEP,FNEXTINFO
C------------------------------------------------------------------------------
        IF(NEXTINFO.GT.10) RETURN
C
        IF(I.EQ.N1)THEN
          WRITE(*,'(A)') 'Working...'
          NEXTINFO=0
        ELSE
          IF(I.LT.N1) RETURN
        END IF
C
        FN1=REAL(N1)
        FN2=REAL(N2)
        FISTEP=REAL(ISTEP)
        FNEXTINFO=REAL(NEXTINFO)
C
        IF((FN2-FN1)/FISTEP.GT.10.)THEN
          IF(I-N1.GE.NINT((FN2-FN1)/10.*FNEXTINFO))THEN
            IF(NEXTINFO.EQ.10)THEN
              WRITE(*,'(A8)'),'100% OK!'
              NEXTINFO=NEXTINFO+1
            ELSE
              WRITE(*,'(I3,A5)')10*NEXTINFO,'%,...'
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
C
        END
