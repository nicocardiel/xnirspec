## Process this file with automake to produce Makefile.in
#
# Unit tests, built and run with "make check". Each test is linked against
# the convenience library ../src/libxnirspec.a (all the routines of xnirspec
# except the main program). A test passes when it exits with status 0.
#
# Run a single test:  make -C tests check TESTS=test_fpoly

check_LIBRARIES = libtestutils.a
libtestutils_a_SOURCES = test_utils.f90

check_PROGRAMS = test_fpoly test_stats test_sort test_linalg test_spline \
                 test_strings readc_file_driver

LDADD = libtestutils.a $(top_builddir)/src/libxnirspec.a

test_fpoly_SOURCES = test_fpoly.f90
test_stats_SOURCES = test_stats.f90
test_sort_SOURCES = test_sort.f90
test_linalg_SOURCES = test_linalg.f90
test_spline_SOURCES = test_spline.f90
test_strings_SOURCES = test_strings.f90
readc_file_driver_SOURCES = readc_file_driver.f90

TESTS = test_fpoly test_stats test_sort test_linalg test_spline \
        test_strings test_readc_file.sh

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)

EXTRA_DIST = test_readc_file.sh

# The tests use the module TEST_UTILS (test_utils.mod), which must be
# compiled first (Automake does not track Fortran module dependencies)
test_fpoly.$(OBJEXT) test_stats.$(OBJEXT) test_sort.$(OBJEXT) \
test_linalg.$(OBJEXT) test_spline.$(OBJEXT) \
test_strings.$(OBJEXT): libtestutils.a

CLEANFILES = *.mod
