C-----------------------------------------------------------------------
C Subrutina extraida de:
C http://light.web.cern.ch/Light/examples/aleph/julia/code/SRCPAT_src.html
      SUBROUTINE SRCPAT(NCOR1,XCOR1,NCOR2,XCOR2,NCOR,XCOM)
C----------------------------------------------------------------------
C! Find intersection between two polygons
C!
C!    Author:     H. Meinhard       26-Feb-1990
C!    Modified:   H. Meinhard       09-Mar-1990  (1)
C!
C!    Input:      - NCOR1     /I    number of corners of 1st polygon
C!                - XCOR1ij   /R    coordinate i of corner j of 1st pl.
C!                - NCOR2     /I    number of corners of 2nd polygon
C!                - XCOR2ij   /R    coordinate i of corner j of 2nd pl.
C!    Output:     - NCOR      /I    number of corners of itersection
C!                - XCOMij    /R    coordinate i of corner j of inters.
C!
C!    Description
C!    ===========
C!    Find the intersection between two polygons being defined by the
C!    number and the coordinates of their corners. The second polygon
C!    must be convex. NCOR2 must be at most 8 (see parameter
C!    statement). The routine uses a clipping algorithm (I. E.
C!    Sutherland, G. W. Hodgman, Communications of the ACM, Vol. 17,
C!    No. 1, January, 1974).
C----------------------------------------------------------------------
*IF .NOT.DOC
      REAL XCOR1(2,*),XCOR2(2,*),XCOM(2,*)
C MXCOR is the maximum number of corners of clipping patch
      PARAMETER (MXCOR = 8)
      REAL F(MXCOR,2),S(MXCOR,2),B1X(MXCOR),B1Y(MXCOR),
     +  B2X(MXCOR),B2Y(MXCOR),B3X(MXCOR),B3Y(MXCOR),ALB3(MXCOR)
      INTEGER IFL(MXCOR),IJUMP(MXCOR)
C----------------------------------------------------------------------
      CALL VZERO(IFL,MXCOR)
      CALL VZERO(IJUMP,MXCOR)
      NCOR = 0
C loop over corners ("vertices") of first polygon
      ICOR1 = 0
  300 ICOR1 = ICOR1 + 1
      IF (ICOR1 .GT. NCOR1)                                 GOTO 350
C loop over clipping straight lines
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
C first point?
      IF (IFL(ILVL) .EQ. 0) THEN
C compute distance measures
        B1X(ILVL) = XCOR2(1,MOD(ILVL-1,NCOR2)+1)
        B1Y(ILVL) = XCOR2(2,MOD(ILVL-1,NCOR2)+1)
        B2X(ILVL) = XCOR2(1,MOD(ILVL,NCOR2)+1)
        B2Y(ILVL) = XCOR2(2,MOD(ILVL,NCOR2)+1)
        B3X(ILVL) = XCOR2(1,MOD(ILVL+1,NCOR2)+1)
        B3Y(ILVL) = XCOR2(2,MOD(ILVL+1,NCOR2)+1)
        ALB3(ILVL) = (B3X(ILVL)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) -
     +               (B3Y(ILVL)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
C P -> S -> F
        IFL(ILVL) = 1
        S(ILVL,1) = PX
        S(ILVL,2) = PY
        F(ILVL,1) = PX
        F(ILVL,2) = PY
      ELSE
C compute distance measures
        ALS = (S(ILVL,1)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) -
     +        (S(ILVL,2)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
        ALP = (PX-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) -
     +        (PY-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
C does line SP cross limiting plane?
        IF (ALS*ALP .LE. 0.) THEN
C compute intersection T of SP and the plane
          IF (ALP-ALS .EQ. 0.) THEN
            TX = PX
            TY = PY
          ELSE
            TX = PX + (ALP/(ALP-ALS))*(S(ILVL,1)-PX)
            TY = PY + (ALP/(ALP-ALS))*(S(ILVL,2)-PY)
          ENDIF
C P -> S
          S(ILVL,1) = PX
          S(ILVL,2) = PY
C output I
          PX = TX
          PY = TY
          IJUMP(ILVL) = 1
          ILVL = ILVL + 1
          GOTO 310
        ELSE
C P -> S
          S(ILVL,1) = PX
          S(ILVL,2) = PY
        ENDIF
      ENDIF
  320   CONTINUE
C is S on visible side of plane?
      ALS = (S(ILVL,1)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) -
     +      (S(ILVL,2)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
      IF (ALS*ALB3(ILVL) .GE. 0.) THEN
        PX = S(ILVL,1)
        PY = S(ILVL,2)
        IJUMP(ILVL) = 2
        ILVL = ILVL + 1
        GOTO 310
      ENDIF
  330   CONTINUE
C exit: go back one level
!!!  340 CONTINUE
      ILVL = ILVL - 1
      IF (ILVL .GE. 1) THEN
        GOTO (320,330,370) IJUMP(ILVL)
      ENDIF
C next vertex
      GOTO 300
C all vertices done
  350 CONTINUE
      ILVL = 1
  360 CONTINUE
C was there any output?
      IF (IJUMP(ILVL) .NE. 0) THEN
        ALS = (S(ILVL,1)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) -
     +        (S(ILVL,2)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
        ALF = (F(ILVL,1)-B1X(ILVL))*(B2Y(ILVL)-B1Y(ILVL)) -
     +        (F(ILVL,2)-B1Y(ILVL))*(B2X(ILVL)-B1X(ILVL))
C does line SF cross limiting plane?
        IF (ALS*ALF .LE. 0.) THEN
C compute intersection T of SF and the plane
          IF (ALF-ALS .EQ. 0.) THEN
            TX = F(ILVL,1)
            TY = F(ILVL,2)
          ELSE
            TX = F(ILVL,1) + (ALF/(ALF-ALS))*(S(ILVL,1)-F(ILVL,1))
            TY = F(ILVL,2) + (ALF/(ALF-ALS))*(S(ILVL,2)-F(ILVL,2))
          ENDIF
C output I
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
C
C******************************************************************************
C imagino que la subrutina VZERO inicializa a cero los elementos de la matriz
C
        SUBROUTINE VZERO(IA,N)
        IMPLICIT NONE
        INTEGER N
        INTEGER IA(N)
C
        INTEGER I
C------------------------------------------------------------------------------
        DO I=1,N
          IA(I)=0
        END DO
        END
