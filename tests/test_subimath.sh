#!/bin/sh
# Test SUBIMATH through subimath_driver, which runs it with buttons in text
# mode and PGPLOT opened with the /NULL device (see the driver source for the
# initial buffers and the expected results of each test case).
#
# Each test case is a list of answers given to SUBIMATH through stdin:
#   @NB   click on button NB (two lines: "y" and the button number)
#   <CR>  empty line
#   other answer to a question (constant, filter option, ...)
#
# Button numbers in SUBIMATH (NORIGEN=86, NMAXBUFF=12):
#   first buffer #k: 86+10*(k-1)     second buffer #k: 89+10*(k-1)
#   operations: 87 [<] copy, 97 [+], 107 [-], 117 [*], 127 [/], 147 [f]ilter
#   region:     88 [f]ull frame, 98 [r]egion
#   operand:    108 cons[t]ant, 118 [f]rame, 128 [x]-cut
#   216 [C]ancel, 217 e[X]it

# "make check" runs the tests from the build directory of tests/, where the
# driver is built; when the script is run by hand from another directory,
# look for the driver next to the script (in-tree builds)
if [ -x ./subimath_driver ]; then
    driver=./subimath_driver
else
    driver=$(dirname "$0")/subimath_driver
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
    out=$(for tok in "$@"; do
              case $tok in
                  (@*)     printf 'y\n%s\n' "${tok#@}" ;;
                  ('<CR>') printf '\n' ;;
                  (*)      printf '%s\n' "$tok" ;;
              esac
          done | "$driver" "$name" 2>&1)
    st=$?
    # If the answers do not match the questions asked by SUBIMATH, the
    # program may end through a STOP statement (exit status 0) before the
    # checks are performed, so the summary line must be present. The prompts
    # are written without a newline, so the results are not always at the
    # beginning of a line.
    if [ $st -eq 0 ] && \
       printf '%s\n' "$out" | grep -q 'Summary: .* 0 failures$'; then
        printf '%s\n' "$out" | grep -oE '(PASS|FAIL): .*' | \
            sed -E "s/^(PASS|FAIL): /\1: [$name] /"
        echo "PASS: case $name"
    else
        printf '%s\n' "$out"
        echo "FAIL: case $name"
        status=1
    fi
}

run_case add_frames          @86 @97  @88 @118 @99 @217
run_case sub_frames          @86 @107 @88 @118 @99 @217
run_case mul_frames          @86 @117 @88 @118 @99 @217
run_case div_frames          @86 @127 @88 @118 @99 @217
# division by zero: [c]ontinue, value -1 for infinity, error -2
run_case div_frames_zero     @86 @127 @88 @118 @99 c -1 -2 @217
# constant 7 with error 0.5
run_case add_constant        @86 @97  @88 @108 7 0.5 @217
# constant 0 is rejected (<CR> to continue; its error is asked anyway),
# then constant 4 with error 0
run_case div_constant_zero   @86 @127 @88 @108 0 '<CR>' 0 4 0 @217
# region X=2..3, Y=2..4, constant 10 with error 0
run_case add_constant_region @86 @97  @98 2 3 2 4 @108 10 0 @217
# buffer #3 <- buffer #2
run_case copy                @106 @87 @99 @217
# filters applied to #1, output in #3
run_case reverse_x           @86 @147 @88 @109 x @217
run_case reverse_y           @86 @147 @88 @109 y @217
run_case median_filter       @86 @147 @88 @109 3 3 3 @217
run_case bin_expand          @86 @147 @88 @109 b x 2 @217
# mean X-cut (1) of #2 with mask value -999, rows 1 to 4 (0 ends the list)
run_case xcut_mask           @86 @107 @88 @128 @99 1 y -999 1 4 0 @217
# as xcut_mask, dividing: column 5 is fully masked (X-cut = 0), so
# [c]ontinue with value -1 for infinity and error -2
run_case xcut_div_zero       @86 @127 @88 @128 @99 1 y -999 1 4 0 c -1 -2 @217

exit $status
