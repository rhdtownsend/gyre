#!/bin/bash
#
# File     : test.sh
# Purpose  : GYRE testing script

. test_support

# Settings

EXEC=./gyre

IN_FILE_N=gyre.n.in
IN_FILE_1=gyre.1.in

OUT_FILE_N=summary.n.h5
OUT_FILE_1=summary.1.h5

# Do the tests

run_gyre $EXEC $IN_FILE_N "numerics (OpenMP)"
if [ $? -ne 0 ]; then
    exit 1;
fi

export OMP_NUM_THREADS=1

run_gyre $EXEC $IN_FILE_1 ""
if [ $? -ne 0 ]; then
    exit 1;
fi

check_output $OUT_FILE_N $OUT_FILE_1
if [ $? -ne 0 ]; then
    exit 1;
fi

# Clean up output files

rm -f $OUT_FILE_N $OUT_FILE_1

# Finish

echo " ...succeeded"
