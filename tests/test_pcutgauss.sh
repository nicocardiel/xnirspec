#!/bin/sh
# Test the interactive gaussian fit PCUTGAUSS through pcutgauss_driver (see
# the driver source for the profile and the expected results). In text
# mode, each click of RPGBAND is given as four lines: "n" (not a button),
# XC, YC and the character (A = left mouse button, X = right button).

# "make check" runs the tests from the build directory of tests/, where the
# driver is built; when the script is run by hand from another directory,
# look for the driver next to the script (in-tree builds)
if [ -x ./pcutgauss_driver ]; then
    driver=./pcutgauss_driver
else
    driver=$(dirname "$0")/pcutgauss_driver
fi
if [ ! -x "$driver" ]; then
    echo "FAIL: $driver not found (run \"make check\" first)"
    exit 1
fi

status=0

# run_case NAME CLICK... (each CLICK is "XC:CH")
run_case() {
    name=$1
    shift
    out=$(for click in "$@"; do
              printf 'n\n%s\n0\n%s\n' "${click%%:*}" "${click##*:}"
          done | "$driver" "$name" 2>&1)
    st=$?
    printf '%s\n' "$out"
    # a mismatched sequence of answers may end the program through STOP
    # (exit status 0) before the checks, so the summary line must be present
    if [ $st -eq 0 ] && \
       printf '%s\n' "$out" | grep -q 'Summary: .* 0 failures$'; then
        echo "PASS: case $name"
    else
        echo "FAIL: case $name"
        status=1
    fi
}

# left continuum 5-25, right continuum 90-55 (any order), line 30-50
run_case fit    5:A 25:A 90:A 55:A 30:A 50:A
# right mouse button at the third click
run_case cancel 5:A 25:A 90:X
# continuum regions without points: the fit must stop after the fourth
# click (if the line region were asked for, stdin would be exhausted)
run_case nocont 200:A 210:A 300:A 310:A

exit $status
