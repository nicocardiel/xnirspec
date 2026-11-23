#!/bin/sh
# Test HPLOT('manual') (button "manual lim.") through hplot_driver (see the
# driver source for the profile and the expected results). Each answer
# (Xmin, Xmax, Ymin, Ymax, in this order) is a line of stdin; an empty line
# keeps the current limit (default value).

# "make check" runs the tests from the build directory of tests/, where the
# driver is built; when the script is run by hand from another directory,
# look for the driver next to the script (in-tree builds)
if [ -x ./hplot_driver ]; then
    driver=./hplot_driver
else
    driver=$(dirname "$0")/hplot_driver
fi
if [ ! -x "$driver" ]; then
    echo "FAIL: $driver not found (run \"make check\" first)"
    exit 1
fi

status=0

# run_case NAME ANSWER... (each ANSWER is a line of stdin; "" = default)
run_case() {
    name=$1
    shift
    if [ $# -gt 0 ]; then
        out=$(printf '%s\n' "$@" | "$driver" "$name" 2>&1)
    else
        out=$("$driver" "$name" 2>&1 < /dev/null)
    fi
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

run_case default   "" "" "" ""
run_case set       10 50 -1 100
run_case retry     50 10 10 50 5 5 0 1
run_case noprofile

exit $status
