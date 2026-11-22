#!/bin/sh
# Test PSAVECUT (save the displayed 1-D profile in a FITS file) through
# psavecut_driver (see the driver source for the test cases). The output
# file name is read from stdin, one answer per line.

# "make check" runs the tests from the build directory of tests/, where the
# driver is built; when the script is run by hand from another directory,
# look for the driver next to the script (in-tree builds)
if [ -x ./psavecut_driver ]; then
    driver=./psavecut_driver
else
    driver=$(dirname "$0")/psavecut_driver
fi
if [ ! -x "$driver" ]; then
    echo "FAIL: $driver not found (run \"make check\" first)"
    exit 1
fi

status=0

# run_case NAME ANSWER...
run_case() {
    name=$1
    shift
    out=$(printf '%s\n' "$@" | "$driver" "$name" 2>&1)
    st=$?
    printf '%s\n' "$out"
    # a mismatched sequence of answers may end the program before the
    # checks, so the summary line must be present
    if [ $st -eq 0 ] && \
       printf '%s\n' "$out" | grep -q 'Summary: .* 0 failures$'; then
        echo "PASS: case $name"
    else
        echo "FAIL: case $name"
        status=1
    fi
}

run_case xwave  tmp_psavecut.fits
run_case ycut   tmp_psavecut.fits
run_case none   none
run_case exists tmp_psavecut_old.fits tmp_psavecut.fits

exit $status
