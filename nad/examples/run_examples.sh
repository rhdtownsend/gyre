#!/bin/sh

for d in \
    mesa/*; do

    cd $d && pwd

    ./gyre_nad < gyre_nad.in | tee gyre_nad.out

    cd -

done
        