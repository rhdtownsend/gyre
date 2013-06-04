#!/bin/sh

for d in \
    mesa/* \
    poly/* \
    hom; do

    cd $d && pwd

    ./gyre_rad gyre_rad.in | tee gyre_rad.out

    cd -

done
        