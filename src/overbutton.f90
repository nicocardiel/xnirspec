!------------------------------------------------------------------------------
! SUBROUTINE OVERBUTTON(LFORCE,NCOLOR)
!
! Draw button 51 ("over=TRUE"/"over=FALSE"), which toggles the overplotting
! of the cuts of other buffers (LOVERCUTS). The button is always visible, but
! it is only active when more than one data buffer (#1 to #NMAXBUFF/2) has
! been defined (NAXIS1 and NAXIS2 different from zero), since otherwise there
! is nothing to overplot. Otherwise it is shown disabled (gray text).
!
! Parameters
! ----------
! LFORCE : LOGICAL
!     If .TRUE., the button is always drawn. If .FALSE., it is only drawn
!     when its state (active/disabled or the value of LOVERCUTS) has changed
!     since the previous call, or in the first call.
! NCOLOR : INTEGER
!     Mode passed to BUTTON for the text of the active button (e.g. NCOLOR4
!     in the main program).
!------------------------------------------------------------------------------
        SUBROUTINE OVERBUTTON(LFORCE,NCOLOR)
        IMPLICIT NONE
        LOGICAL LFORCE
        INTEGER NCOLOR
!
        INCLUDE 'dimensions.inc'
!
        INTEGER I,NLOADED
        INTEGER NAXIS(2,NMAXBUFF)
        CHARACTER*10 CLABEL
        LOGICAL LACTIVE
        LOGICAL LOVERCUTS
        LOGICAL, SAVE :: LFIRST=.TRUE.
        LOGICAL, SAVE :: LACTIVE_OLD=.FALSE.
        LOGICAL, SAVE :: LOVERCUTS_OLD=.FALSE.
!
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKLOVERCUTS/LOVERCUTS
!------------------------------------------------------------------------------
! number of defined data buffers (same criterion as in SUBIMATH)
        NLOADED=0
        DO I=1,NMAXBUFF/2
          IF((NAXIS(1,I).NE.0).AND.(NAXIS(2,I).NE.0)) NLOADED=NLOADED+1
        END DO
        LACTIVE=(NLOADED.GT.1)
!
        IF(LFORCE.OR.LFIRST.OR.(LACTIVE.NEQV.LACTIVE_OLD).OR. &
         (LOVERCUTS.NEQV.LOVERCUTS_OLD))THEN
          IF(LOVERCUTS)THEN
            CLABEL='over=TRUE'
          ELSE
            CLABEL='over=FALSE'
          END IF
          CALL BUTTON(51,TRIM(CLABEL),0)
          IF(LACTIVE)THEN
            CALL BUTTON(51,TRIM(CLABEL),NCOLOR)
          ELSE
            CALL BUTTON(51,TRIM(CLABEL),3)                 !disabled button
          END IF
          LFIRST=.FALSE.
          LACTIVE_OLD=LACTIVE
          LOVERCUTS_OLD=LOVERCUTS
        END IF
        END
