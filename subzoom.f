        SUBROUTINE SUBZOOM(NX1,NX2,NY1,NY2,LOK)
        IMPLICIT NONE
        INTEGER NX1,NX2,NY1,NY2
        LOGICAL LOK
C
        INTEGER IXC1,IXC2,IYC1,IYC2
        REAL XC,YC
        CHARACTER*1 CH
C------------------------------------------------------------------------------
        LOK=.TRUE.
        WRITE(*,101) 'Press cursor at two corners of the imaginary '
     +   //'BOX to be zoomed'
        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        IF(CH.EQ.'X')THEN
          LOK=.FALSE.
          RETURN
        END IF
        IXC1=INT(XC+0.5)
        IF(IXC1.LT.NX1) IXC1=NX1
        IF(IXC1.GT.NX2) IXC1=NX2
        IYC1=INT(YC+0.5)
        IF(IYC1.LT.NY1) IYC1=NY1
        IF(IYC1.GT.NY2) IYC1=NY2
        WRITE(*,112) 'Cursor at ',IXC1,IYC1
C
        CALL PGSCI(5)
        CALL RPGBAND(2,0,REAL(IXC1),REAL(IYC1),XC,YC,CH)
        CALL PGSCI(1)
        IF(CH.EQ.'X')THEN
          LOK=.FALSE.
          RETURN
        END IF
        IXC2=INT(XC+0.5)
        IF(IXC2.LT.NX1) IXC2=NX1
        IF(IXC2.GT.NX2) IXC2=NX2
        IYC2=INT(YC+0.5)
        IF(IYC2.LT.NY1) IYC2=NY1
        IF(IYC2.GT.NY2) IYC2=NY2
        WRITE(*,112) 'Cursor at ',IXC2,IYC2
C
        NX1=MIN0(IXC1,IXC2)
        NX2=MAX0(IXC1,IXC2)
        NY1=MIN0(IYC1,IYC2)
        NY2=MAX0(IYC1,IYC2)
C
101     FORMAT(A)
112     FORMAT(A,I5,2X,I5)
        END
