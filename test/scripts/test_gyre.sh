#!/bin/bash
#
# File     : test_gyre.sh
# Purpose  : GYRE testing script

# Arguments

EXEC=$1
FILE=$2

# Run GYRE

$EXEC $FILE

if [ $? -ne 0 ]; then
    echo "Failed to run $EXEC"
    exit 1
fi

# Check the output

for ref_file in ref/*.txt; do

    file=${ref_file##ref/}

    ndiff $file $ref_file

    if [ $? -ne 0 ]; then
	echo "Failed when comparing $file and $ref_file"
	exit 1
    fi

done

# Clean up output files

rm -f *.txt

# Finish

echo "Tests passed"
