#!/bin/sh
# Test the paths of SLEEFITS and SESCRFITS that read answers from stdin
# (3D cubes, files with extensions, header templates) through
# fitsio_interactive_driver, which performs the checks itself. See the
# driver source for the meaning of each answer.

# "make check" runs the tests from the build directory of tests/, where the
# driver is built; when the script is run by hand from another directory,
# look for the driver next to the script (in-tree builds)
if [ -x ./fitsio_interactive_driver ]; then
    driver=./fitsio_interactive_driver
else
    driver=$(dirname "$0")/fitsio_interactive_driver
fi
if [ ! -x "$driver" ]; then
    echo "FAIL: $driver not found (run \"make check\" first)"
    exit 1
fi

printf '%s\n' 0 2 -1 2 3 2 '' '' '' | "$driver"
