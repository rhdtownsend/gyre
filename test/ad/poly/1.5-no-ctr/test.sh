#!/bin/bash
#
# File     : test.sh
# Purpose  : GYRE testing script

. test_support

# Settings

EXEC=./gyre

IN_FILE=gyre.in
OUT_FILE=summary.txt

LABEL="polytrope model (n_poly=1.5, no center)"

ABSERR=1E-6
FIELDS=1-2,5

# Do the tests

run_gyre $EXEC $IN_FILE "$LABEL"
if [ $? -ne 0 ]; then
    exit 1;
fi

check_output $ABSERR $FIELDS $OUT_FILE '' abs
if [ $? -ne 0 ]; then
    exit 1;
fi

# Clean up output files

rm -f $OUT_FILE

# Finish

echo " ...succeeded"
