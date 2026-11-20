#!/bin/sh
# Test READC_FILE (fallback path without readline, since stdin is not a
# terminal) through readc_file_driver.
#
# Input lines: empty (-> default value), an explicit file name, an answer
# with invalid characters (rejected) and a valid one.

# "make check" runs the tests from the build directory of tests/, where the
# driver is built; when the script is run by hand from another directory,
# look for the driver next to the script (in-tree builds)
if [ -x ./readc_file_driver ]; then
    driver=./readc_file_driver
else
    driver=$(dirname "$0")/readc_file_driver
fi
if [ ! -x "$driver" ]; then
    echo "FAIL: $driver not found (run \"make check\" first)"
    exit 1
fi

out=$(printf '\n%s\n%s\n%s\n' 'data/my_file.fits' 'abx' 'cab' | \
      "$driver") || { echo "FAIL: readc_file_driver exited with error"; exit 1; }

expected='RESULT:default.fits|
RESULT:data/my_file.fits|
RESULT:cab|'

# the prompts are written without a newline, so RESULT: is not at the
# beginning of the line
got=$(printf '%s\n' "$out" | grep -o 'RESULT:[^|]*|')
nerr=$(printf '%s\n' "$out" | grep -c 'ERROR: invalid character')

status=0
if [ "$got" = "$expected" ]; then
    echo "PASS: values returned by READC_FILE"
else
    echo "FAIL: values returned by READC_FILE"
    echo "--- got:"; printf '%s\n' "$got"
    echo "--- expected:"; printf '%s\n' "$expected"
    status=1
fi
if [ "$nerr" -eq 1 ]; then
    echo "PASS: invalid answer rejected once"
else
    echo "FAIL: invalid answer rejected $nerr times (expected 1)"
    status=1
fi
exit $status
