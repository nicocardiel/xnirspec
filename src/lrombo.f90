! comprueba si un punto XC,YC estan dentro del rombo definido por los
! cuatro vertices introducidos
        LOGICAL FUNCTION LROMBO(XC,YC,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
        IMPLICIT NONE
        REAL XC,YC
        REAL X1,Y1
        REAL X2,Y2
        REAL X3,Y3
        REAL X4,Y4
!
        INTEGER CONTROL
        REAL VPX,VPY
        REAL VVX,VVY
!------------------------------------------------------------------------------
        CONTROL=0
!..............................................................................
        VPX=XC-X1
        VPY=YC-Y1
        VVX=X2-X1
        VVY=Y2-Y1
        IF((VPX*VVY-VVX*VPY).LE.0.) CONTROL=CONTROL+1
!..............................................................................
        VPX=XC-X2
        VPY=YC-Y2
        VVX=X3-X2
        VVY=Y3-Y2
        IF((VPX*VVY-VVX*VPY).LE.0.) CONTROL=CONTROL+10
!..............................................................................
        VPX=XC-X3
        VPY=YC-Y3
        VVX=X4-X3
        VVY=Y4-Y3
        IF((VPX*VVY-VVX*VPY).LE.0.) CONTROL=CONTROL+100
!..............................................................................
        VPX=XC-X4
        VPY=YC-Y4
        VVX=X1-X4
        VVY=Y1-Y4
        IF((VPX*VVY-VVX*VPY).LE.0.) CONTROL=CONTROL+1000
!..............................................................................
        LROMBO=((CONTROL.EQ.0000).OR.(CONTROL.EQ.1111))
!
        END
