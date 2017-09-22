#!/bin/bash
#
# File     : test.sh
# Purpose  : GYRE testing script

. test_support

# Settings

EXEC=./gyre

IN_FILE=gyre.in
OUT_FILE=summary.txt

LABEL="MESA model for slowly pulsating B-type star (adiabatic, quasi-adiabatic eigfuncs)"

RELERR=1E-13
FIELDS=1-7

# Do the tests

run_gyre $EXEC $IN_FILE "$LABEL"
if [ $? -ne 0 ]; then
    exit 1;
fi

check_output $RELERR $FIELDS $OUT_FILE 
if [ $? -ne 0 ]; then
    exit 1;
fi

# Clean up output files

rm -f $OUT_FILE

# Finish

echo " ...succeeded"
