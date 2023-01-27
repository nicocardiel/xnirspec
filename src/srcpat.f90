!-----------------------------------------------------------------------
! Subrutina extraida de:
! http://light.web.cern.ch/Light/examples/aleph/julia/code/SRCPAT_src.html
      SUBROUTINE SRCPAT(NCOR1,XCOR1,NCOR2,XCOR2,NCOR,XCOM)
!----------------------------------------------------------------------
!! Find intersection between two polygons
!!
!!    Author:     H. Meinhard       26-Feb-1990
!!    Modified:   H. Meinhard       09-Mar-1990  (1)
!!
!!    Input:      - NCOR1     /I    number of corners of 1st polygon
!!                - XCOR1ij   /R    coordinate i of corner j of 1st pl.
!!                - NCOR2     /I    number of corners of 2nd polygon
!!                - XCOR2ij   /R    coordinate i of corner j of 2nd pl.
!!    Output:     - NCOR      /I    number of corners of itersection
!!                - XCOMij    /R    coordinate i of corner j of inters.
!!
!!    Description
!!    ===========
!!    Find the intersection between two polygons being defined by the
!!    number and the coordinates of their corners. The second polygon
!!    must be convex. NCOR2 must be at most 8 (see parameter
!!    statement). The routine uses a clipping algorithm (I. E.
!!    Sutherland, G. W. Hodgman, Communications of the ACM, Vol. 17,
!!    No. 1, January, 1974).
!----------------------------------------------------------------------
!*IF .NOT.DOC
      REAL XCOR1(2,*),XCOR2(2,*),XCOM(2,*)
! MXCOR is the maximum number of corners of clipping patch
      PARAMETER (MXCOR = 8)
      REAL F(MXCOR,2),S(MXCOR,2),B1X(MXCOR),B1Y(MXCOR),B2X(MXCOR),B2Y(MXCOR),B3X(MXCOR),B3Y(MXCOR),ALB3(MXCOR)
      INTEGER IFL(MXCOR),IJUMP(MXCOR)
!----------------------------------------------------------------------
      CALL VZERO(IFL,MXCOR)
      CALL VZERO(IJUMP,MXCOR)
      NCOR = 0
! loop over corners ("vertices") of first polygon
      ICOR1 = 0
  300 ICOR1 = ICOR1 + 1
      IF (ICOR1 .GT. NCOR1)                                 GOTO 350
! loop over clipping straight lines
      ILVL = 1
      PX = XCOR1(1,ICOR1)
      PY = XCOR1(2,ICOR1)
  310 CONTINUE
      IF (ILVL .EQ. NCOR2 + 1) THEN
        NCOR = NCOR + 1
        XCOM(1,NCOR) = PX
        XCOM(2,NCOR) = PY
        GOTO 330
      ENDIF
! first point?
      IF (IFL(ILVL) .EQ. 0) THEN
! compute distance measures
        B1X(ILVL) = XCOR2(1,MOD(ILVL-1,NCOR2)+1)
        B1Y(ILVL) = XCOR2(2,MOD(ILVL-1,NCOR2)+1)
        B2X(ILVL) = XCOR2(1,MOD(ILVL,NCOR2)+1)
        B2Y(ILVL) = XCOR2(2,MOD(ILVL,NCOR2)+1)
        B3X(ILVL) = XCOR2(1,MOD(ILVL+1,NCOR2)+1)
        B3Y(ILVL) = XCOR2(2,MOD(ILVL+1,NCOR2)+1)
        ALB3(ILVL) = (B3X(ILVL)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) - (B3Y(ILVL)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
! P -> S -> F
        IFL(ILVL) = 1
        S(ILVL,1) = PX
        S(ILVL,2) = PY
        F(ILVL,1) = PX
        F(ILVL,2) = PY
      ELSE
! compute distance measures
        ALS = (S(ILVL,1)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) - (S(ILVL,2)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
        ALP = (PX-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) - (PY-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
! does line SP cross limiting plane?
        IF (ALS*ALP .LE. 0.) THEN
! compute intersection T of SP and the plane
          IF (ALP-ALS .EQ. 0.) THEN
            TX = PX
            TY = PY
          ELSE
            TX = PX + (ALP/(ALP-ALS))*(S(ILVL,1)-PX)
            TY = PY + (ALP/(ALP-ALS))*(S(ILVL,2)-PY)
          ENDIF
! P -> S
          S(ILVL,1) = PX
          S(ILVL,2) = PY
! output I
          PX = TX
          PY = TY
          IJUMP(ILVL) = 1
          ILVL = ILVL + 1
          GOTO 310
        ELSE
! P -> S
          S(ILVL,1) = PX
          S(ILVL,2) = PY
        ENDIF
      ENDIF
  320   CONTINUE
! is S on visible side of plane?
      ALS = (S(ILVL,1)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) - (S(ILVL,2)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
      IF (ALS*ALB3(ILVL) .GE. 0.) THEN
        PX = S(ILVL,1)
        PY = S(ILVL,2)
        IJUMP(ILVL) = 2
        ILVL = ILVL + 1
        GOTO 310
      ENDIF
  330   CONTINUE
! exit: go back one level
!!!  340 CONTINUE
      ILVL = ILVL - 1
      IF (ILVL .GE. 1) THEN
        GOTO (320,330,370) IJUMP(ILVL)
      ENDIF
! next vertex
      GOTO 300
! all vertices done
  350 CONTINUE
      ILVL = 1
  360 CONTINUE
! was there any output?
      IF (IJUMP(ILVL) .NE. 0) THEN
        ALS = (S(ILVL,1)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) - (S(ILVL,2)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
        ALF = (F(ILVL,1)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) - (F(ILVL,2)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
! does line SF cross limiting plane?
        IF (ALS*ALF .LE. 0.) THEN
! compute intersection T of SF and the plane
          IF (ALF-ALS .EQ. 0.) THEN
            TX = F(ILVL,1)
            TY = F(ILVL,2)
          ELSE
            TX = F(ILVL,1) + (ALF/(ALF-ALS))*(S(ILVL,1)-F(ILVL,1))
            TY = F(ILVL,2) + (ALF/(ALF-ALS))*(S(ILVL,2)-F(ILVL,2))
          ENDIF
! output I
          PX = TX
          PY = TY
          IJUMP(ILVL) = 3
          ILVL = ILVL + 1
          GOTO 310
        ENDIF
      ENDIF
  370   CONTINUE
      IFL(ILVL) = 0
      IF (ILVL .LT. NCOR2) THEN
        ILVL = ILVL + 1
        GOTO 360
      ENDIF
      END
!
!******************************************************************************
! imagino que la subrutina VZERO inicializa a cero los elementos de la matriz
!
        SUBROUTINE VZERO(IA,N)
        IMPLICIT NONE
        INTEGER N
        INTEGER IA(N)
!
        INTEGER I
!------------------------------------------------------------------------------
        DO I=1,N
          IA(I)=0
        END DO
        END
