!------------------------------------------------------------------------------
! Unit tests for TRUELEN, TRUEBEG, RMBLANK and CHUPPER.
!------------------------------------------------------------------------------
        PROGRAM TEST_STRINGS
        USE TEST_UTILS
        IMPLICIT NONE
!
        INTEGER TRUELEN,TRUEBEG
        INTEGER L
        CHARACTER*20 C20
        CHARACTER*9 C9
!------------------------------------------------------------------------------
        CALL CHECK_INT('TRUELEN leading/trailing blanks',TRUELEN('  abc  '),5)
        CALL CHECK_INT('TRUELEN blank string',TRUELEN('   '),0)
        CALL CHECK_INT('TRUELEN zero-length string',TRUELEN(''),0)
        CALL CHECK_INT('TRUELEN no blanks',TRUELEN('abc'),3)
!
        CALL CHECK_INT('TRUEBEG leading/trailing blanks',TRUEBEG('  abc  '),3)
        CALL CHECK_INT('TRUEBEG blank string',TRUEBEG('   '),0)
        CALL CHECK_INT('TRUEBEG zero-length string',TRUEBEG(''),0)
        CALL CHECK_INT('TRUEBEG no blanks',TRUEBEG('abc'),1)
!------------------------------------------------------------------------------
        C20='xxxxxxxxxxxxxxxxxxxx'
        CALL RMBLANK(' a b  c ',C20,L)
        CALL CHECK_STR('RMBLANK string',C20,'abc')
        CALL CHECK_INT('RMBLANK length',L,3)
        CALL RMBLANK('    ',C20,L)
        CALL CHECK_INT('RMBLANK blank string length',L,0)
!------------------------------------------------------------------------------
        C9='abc-Xyz 1'
        CALL CHUPPER(C9)
        CALL CHECK_STR('CHUPPER',C9,'ABC-XYZ 1')
!
        CALL TEST_SUMMARY()
        END PROGRAM TEST_STRINGS
