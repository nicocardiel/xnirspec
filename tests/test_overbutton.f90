!------------------------------------------------------------------------------
! Unit tests for OVERBUTTON: button 51 (overplot of cuts) must be active
! (EXIST_BUTT(51)=.TRUE.) only when more than one data buffer (#1 to
! #NMAXBUFF/2) is defined; button 44 (gaussian fit to the displayed cut)
! and button 54 (save the displayed cut) only when at least one data
! buffer is defined and LOVERCUTS is .FALSE.
! Buttons are used in text mode without plotting them, so no graphic device
! is needed.
!------------------------------------------------------------------------------
        PROGRAM TEST_OVERBUTTON
        USE TEST_UTILS
        IMPLICIT NONE
        INCLUDE 'dimensions.inc'
        INCLUDE 'button.inc'
!
        INTEGER NCOLOR4
        PARAMETER (NCOLOR4=-5-1)                     !as in xnirspec.f90
! global variables (COMMONs) shared with OVERBUTTON
        INTEGER NAXIS(2,NMAXBUFF)
        LOGICAL LOVERCUTS
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKLOVERCUTS/LOVERCUTS
!------------------------------------------------------------------------------
        MODOTEXT_BUTT=.TRUE.
        MODOTEXT_PLOTBUTT=.FALSE.
        EXIST_BUTT=.FALSE.
        CALL BUTTSXB(10)
        CALL BUTTSYB(25)
        LOVERCUTS=.FALSE.
        NAXIS=0
! no buffer: disabled (the first call always draws the button)
        CALL OVERBUTTON(.FALSE.,NCOLOR4)
        CALL CHECK_TRUE('no buffer defined: disabled',.NOT.EXIST_BUTT(51))
        CALL CHECK_TRUE('no buffer defined: 44 disabled',.NOT.EXIST_BUTT(44))
        CALL CHECK_TRUE('no buffer defined: 54 disabled',.NOT.EXIST_BUTT(54))
! one buffer: disabled
        NAXIS(:,1)=(/5,4/)
        CALL OVERBUTTON(.FALSE.,NCOLOR4)
        CALL CHECK_TRUE('one buffer defined: disabled',.NOT.EXIST_BUTT(51))
        CALL CHECK_TRUE('one buffer defined: 44 active',EXIST_BUTT(44))
        CALL CHECK_TRUE('one buffer defined: 54 active',EXIST_BUTT(54))
! two buffers: active (redrawn because the state changed)
        NAXIS(:,3)=(/7,2/)
        CALL OVERBUTTON(.FALSE.,NCOLOR4)
        CALL CHECK_TRUE('two buffers defined: active',EXIST_BUTT(51))
! no change and LFORCE=.FALSE.: the button is not redrawn
        EXIST_BUTT(51)=.FALSE.
        CALL OVERBUTTON(.FALSE.,NCOLOR4)
        CALL CHECK_TRUE('no change: not redrawn',.NOT.EXIST_BUTT(51))
! LFORCE=.TRUE.: always redrawn
        CALL OVERBUTTON(.TRUE.,NCOLOR4)
        CALL CHECK_TRUE('LFORCE: redrawn (active)',EXIST_BUTT(51))
! toggling LOVERCUTS with two buffers: still active
        LOVERCUTS=.TRUE.
        CALL OVERBUTTON(.FALSE.,NCOLOR4)
        CALL CHECK_TRUE('LOVERCUTS toggled: active',EXIST_BUTT(51))
        CALL CHECK_TRUE('LOVERCUTS=.TRUE.: 44 disabled',.NOT.EXIST_BUTT(44))
        CALL CHECK_TRUE('LOVERCUTS=.TRUE.: 54 disabled',.NOT.EXIST_BUTT(54))
! error buffers do not count: data buffer #1 and error buffer #3
        NAXIS(:,3)=0
        NAXIS(:,3+NMAXBUFF/2)=(/7,2/)
        CALL OVERBUTTON(.FALSE.,NCOLOR4)
        CALL CHECK_TRUE('error buffers not counted: disabled', &
         .NOT.EXIST_BUTT(51))
! back to LOVERCUTS=.FALSE. with one data buffer: 44 active again
        LOVERCUTS=.FALSE.
        CALL OVERBUTTON(.FALSE.,NCOLOR4)
        CALL CHECK_TRUE('LOVERCUTS=.FALSE.: 44 active',EXIST_BUTT(44))
        CALL CHECK_TRUE('LOVERCUTS=.FALSE.: 54 active',EXIST_BUTT(54))
!
        CALL TEST_SUMMARY()
        END PROGRAM TEST_OVERBUTTON
