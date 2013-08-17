#!/bin/bash
#
# File     : test_gyre.sh
# Purpose  : GYRE testing script

# Arguments

EXEC=$1
FILE=$2
LABEL=$3

# Run GYRE

echo -n "TEST $LABEL"

$EXEC $FILE > /dev/null

if [ $? -ne 0 ]; then
    echo " ...failed during execution of $EXEC"
    exit 1
fi

# Check the output

for ref_file in ref/*.txt; do

    file=${ref_file##ref/}

    ndiff -quiet -relerr 1E-12 -fields 1-5 $file $ref_file

    if [ $? -ne 0 ]; then
	echo " ...failed when comparing $file and $ref_file"
	exit 1
    fi

done

# Clean up output files

rm -f *.txt

# Finish

echo " ...succeeded"
