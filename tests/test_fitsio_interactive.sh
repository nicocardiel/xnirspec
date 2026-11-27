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

output=$(printf '%s\n' 0 2 -1 2 3 -1 2 0 1 '' '' '' '' '' | "$driver")
status=$?
printf '%s\n' "$output"
if [ $status -ne 0 ]; then
    exit $status
fi

# a mismatched answer sequence can end the driver through STOP (exit status 0)
if ! printf '%s\n' "$output" | grep -q 'Summary: .* 0 failures'; then
    echo "FAIL: the driver did not finish with 0 failures"
    exit 1
fi

# listing of the HDUs printed with the answer -1 (astropy fitsinfo layout)
for pattern in \
    '^No\.    Name      Ver    Type      Cards   Dimensions   Format$' \
    '^  0  PRIMARY       1 PrimaryHDU  ' \
    '^  1  SCI           1 ImageHDU    .*   (4, 3)   float32$' \
    '^  2  TAB           1 BinTableHDU .*   5R x 2C   \[E, J\]$'; do
    if ! printf '%s\n' "$output" | grep -q "$pattern"; then
        echo "FAIL: HDU listing line not found: $pattern"
        exit 1
    fi
done
