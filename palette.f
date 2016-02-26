        SUBROUTINE PALETTE(ITYPE)
        IMPLICIT NONE
        INTEGER ITYPE
C
        REAL CONTRA, BRIGHT
C
        REAL GL(2), GR(2), GG(2), GB(2)
        REAL RL(9), RR(9), RG(9), RB(9)
        REAL HL(5), HR(5), HG(5), HB(5)
        REAL WL(10), WR(10), WG(10), WB(10)
        REAL AL(20), AR(20), AG(20), AB(20)
C
        DATA GL /0.0, 1.0/
        DATA GR /0.0, 1.0/
        DATA GG /0.0, 1.0/
        DATA GB /0.0, 1.0/
C
        DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
        DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
        DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
        DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
        DATA HL /0.0, 0.3, 0.6, 0.9, 1.0/
        DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
        DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
        DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
        DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
        DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
        DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
        DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
        DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     +           0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
        DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     +           0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
        DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     +           0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
        DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     +           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C------------------------------------------------------------------------------
        CONTRA=1.0
        CONTRA=-CONTRA
        BRIGHT=0.5
        BRIGHT=1-BRIGHT
C
        IF (ITYPE.EQ.1) THEN
C         -- gray scale
          CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
          WRITE(*,101) '>>> palette: gray scale'
        ELSE IF (ITYPE.EQ.2) THEN
C         -- rainbow
          CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
          WRITE(*,101) '>>> palette: rainbow'
        ELSE IF (ITYPE.EQ.3) THEN
C         -- heat
          CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
          WRITE(*,101) '>>> palette: heat'
        ELSE IF (ITYPE.EQ.4) THEN
C         -- weird IRAF
          CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
          WRITE(*,101) '>>> palette: weird IRAF'
        ELSE IF (ITYPE.EQ.5) THEN
C         -- AIPS
          CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
          WRITE(*,101) '>>> palette: AIPS'
        END IF
C------------------------------------------------------------------------------
101     FORMAT(A)
        END
