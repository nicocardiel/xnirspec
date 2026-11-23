# Process this file with autoconf to produce a configure script.
#

# Prelude
AC_INIT(xnirspec,6.6,[cardiel@ucm.es])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs
AC_PROG_FC
AC_PROG_CC
AM_PROG_AR
AC_PROG_RANLIB

# Safety check
AC_CONFIG_SRCDIR([src/xnirspec.f90])

# Define directories to search for libraries
AC_MSG_CHECKING([uname -s for detecting host operating system])
OSNAME=`uname -s`
case "$OSNAME" in
  Darwin*)
    AC_MSG_RESULT(Darwin)
    MY_DIR="/opt/local/lib"
    if test -d "$MY_DIR"; then
      LDFLAGS="$LDFLAGS -L$MY_DIR"
    fi
  ;;
  Linux*)
    AC_MSG_RESULT(Linux)
    MY_DIR="/usr/local/pgplot"
    if test -d "$MY_DIR"; then
      LDFLAGS="$LDFLAGS -L$MY_DIR"
    fi
  ;;
esac

# Checks for libraries
AC_MSG_RESULT(LDFLAGS set to... $LDFLAGS)
AC_CHECK_LIB([X11], main,,AC_MSG_ERROR(X11 library not found))
AC_CHECK_LIB([pgplot], main,,AC_MSG_ERROR(pgplot library not found))

#PKG_CHECK_MODULES([CFITSIO], [cfitsio])
AC_CHECK_LIB([cfitsio], main,,AC_MSG_ERROR(cfitsio library not found))

# Optional GNU readline (file name completion in READC_FILE).
# --with-readline=DIR looks for DIR/include and DIR/lib.
AC_ARG_WITH([readline],
  [AS_HELP_STRING([--with-readline@<:@=DIR@:>@],
    [use readline for file name completion @<:@default=check@:>@])],
  [], [with_readline=check])
have_readline=no
AS_IF([test "x$with_readline" != xno], [
  rl_save_CPPFLAGS=$CPPFLAGS
  rl_save_LDFLAGS=$LDFLAGS
  rl_save_LIBS=$LIBS
  case "$with_readline" in
    yes|check)
      # Homebrew installs GNU readline as keg-only (not linked into the
      # default search paths, where the libedit headers of the SDK are found);
      # MacPorts installs it under /opt/local
      case "$OSNAME" in
        Darwin*)
          for rl_dir in /opt/homebrew/opt/readline /usr/local/opt/readline \
                        /opt/local; do
            if test -f "$rl_dir/include/readline/readline.h"; then
              CPPFLAGS="$CPPFLAGS -I$rl_dir/include"
              LDFLAGS="$LDFLAGS -L$rl_dir/lib"
              break
            fi
          done
        ;;
      esac
    ;;
    *)
      CPPFLAGS="$CPPFLAGS -I$with_readline/include"
      LDFLAGS="$LDFLAGS -L$with_readline/lib"
    ;;
  esac
  AC_LANG_PUSH([C])
  AC_CHECK_HEADER([readline/readline.h], [rl_header=yes], [rl_header=no],
    [#include <stdio.h>])
  AS_IF([test "x$rl_header" = xyes], [
    # some readline builds are not linked against a termcap library
    for rl_extra in "" -ltinfo -lncurses -ltermcap; do
      AC_MSG_CHECKING([for readline with "$rl_extra"])
      LIBS="-lreadline $rl_extra $rl_save_LIBS"
      AC_LINK_IFELSE([AC_LANG_CALL([], [readline])],
        [have_readline=yes])
      AC_MSG_RESULT([$have_readline])
      test "x$have_readline" = xyes && break
    done
  ])
  AS_IF([test "x$have_readline" = xyes], [
    AC_DEFINE([HAVE_READLINE], [1],
      [Define to 1 if readline is available.])
    AC_CHECK_HEADERS([readline/history.h], [], [], [#include <stdio.h>])
  ], [
    CPPFLAGS=$rl_save_CPPFLAGS
    LDFLAGS=$rl_save_LDFLAGS
    LIBS=$rl_save_LIBS
    AS_IF([test "x$with_readline" != xcheck],
      [AC_MSG_ERROR([readline requested but not found])])
  ])
  AC_LANG_POP([C])
])
AC_MSG_NOTICE([file name completion with readline: $have_readline])

# Define in LIBS all the previous libraries
LIBS="$LIBS"
AC_MSG_RESULT(LIBS set to... $LIBS)

# Output files
AC_CONFIG_FILES([Makefile 
                 src/Makefile
                 tests/Makefile
                 src/dimensions.inc
                 src/largest.inc
                 src/sexpath.inc
                 src/version.inc])

# Definitions
DEFAULT_NXMAX="4300"
NXMAX=${NXMAX="$DEFAULT_NXMAX"}
AC_ARG_VAR(NXMAX,[Maximum image dimension in the X direction])

DEFAULT_NYMAX="4300"
NYMAX=${NYMAX="$DEFAULT_NYMAX"}
AC_ARG_VAR(NYMAX,[Maximum image dimension in the Y direction])

DEFAULT_NXYMAX="4300"
NXYMAX=${NXYMAX="$DEFAULT_NXYMAX"}
AC_ARG_VAR(NXYMAX,[Maximum image dimension in X or Y direction])

DEFAULT_SEXPATH="/opt/local/bin/sex"
SEXPATH=${SEXPATH="$DEFAULT_SEXPATH"}
AC_ARG_VAR(SEXPATH,[Sextractor path])

DEFAULT_NXMAXB9="900"
NXMAXB9=${NXMAXB9="$DEFAULT_NXMAXB9"}
AC_ARG_VAR(NXMAXB9,[nxmaxb9 value (do not change it!)])

DEFAULT_NYMAXB9="900"
NYMAXB9=${NYMAXB9="$DEFAULT_NYMAXB9"}
AC_ARG_VAR(NYMAXB9,[nymaxb9 value (do not change it!)])

DEFAULT_NOVERSAMPMAX="10"
NOVERSAMPMAX=${NOVERSAMPMAX="$DEFAULT_NOVERSAMPMAX"}
AC_ARG_VAR(NOVERSAMPMAX,[noversampmax value (do not change it!)])

DEFAULT_NZMAX="10000"
NZMAX=${NZMAX="$DEFAULT_NZMAX"}
AC_ARG_VAR(NZMAX,[nzmax value (do not change it!)])

DEFAULT_NMAXBUFF="12"
NMAXBUFF=${NMAXBUFF="$DEFAULT_NMAXBUFF"}
AC_ARG_VAR(NMAXBUFF,[nmaxbuff value (do not change it!)])

AC_OUTPUT

