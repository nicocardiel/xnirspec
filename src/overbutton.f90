!------------------------------------------------------------------------------
! SUBROUTINE OVERBUTTON(LFORCE,NCOLOR)
!
! Draw the buttons that depend on the cuts displayed in the upper-left plot:
!
! - Button 51 ("over=TRUE"/"over=FALSE") toggles the overplotting of the
!   cuts of other buffers (LOVERCUTS). It is always visible, but it is only
!   active when more than one data buffer (#1 to #NMAXBUFF/2) has been
!   defined (NAXIS1 and NAXIS2 different from zero), since otherwise there
!   is nothing to overplot.
! - Button 44 ("gauss fit") fits a gaussian to the displayed 1-D profile
!   (PCUTGAUSS). It is always visible, but it is only active when at least
!   one data buffer has been defined and LOVERCUTS is .FALSE. (i.e. a single
!   profile is displayed).
! - Button 54 ("save cut") saves the displayed 1-D profile in a FITS file
!   (PSAVECUT). It is active under the same conditions as button 44.
!
! Disabled buttons are shown with gray text.
!
! Parameters
! ----------
! LFORCE : LOGICAL
!     If .TRUE., the buttons are always drawn. If .FALSE., they are only
!     drawn when their state (active/disabled or the value of LOVERCUTS) has
!     changed since the previous call, or in the first call.
! NCOLOR : INTEGER
!     Mode passed to BUTTON for the text of the active buttons (e.g.
!     NCOLOR4 in the main program).
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
        LOGICAL LACTIVE,LFITACTIVE
        LOGICAL LOVERCUTS
        LOGICAL, SAVE :: LFIRST=.TRUE.
        LOGICAL, SAVE :: LACTIVE_OLD=.FALSE.
        LOGICAL, SAVE :: LFITACTIVE_OLD=.FALSE.
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
        LFITACTIVE=((NLOADED.GE.1).AND.(.NOT.LOVERCUTS))
!
        IF(LFORCE.OR.LFIRST.OR.(LACTIVE.NEQV.LACTIVE_OLD).OR. &
         (LFITACTIVE.NEQV.LFITACTIVE_OLD).OR. &
         (LOVERCUTS.NEQV.LOVERCUTS_OLD))THEN
! button 51: overplot of cuts
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
! button 44: gaussian fit to the displayed profile
          CALL BUTTON(44,'gauss fit',0)
          IF(LFITACTIVE)THEN
            CALL BUTTON(44,'gauss fit',NCOLOR)
          ELSE
            CALL BUTTON(44,'gauss fit',3)                  !disabled button
          END IF
! button 54: save the displayed profile in a FITS file
          CALL BUTTON(54,'save cut',0)
          IF(LFITACTIVE)THEN
            CALL BUTTON(54,'save cut',NCOLOR)
          ELSE
            CALL BUTTON(54,'save cut',3)                   !disabled button
          END IF
          LFIRST=.FALSE.
          LACTIVE_OLD=LACTIVE
          LFITACTIVE_OLD=LFITACTIVE
          LOVERCUTS_OLD=LOVERCUTS
        END IF
        END
