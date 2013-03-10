#!/bin/sh

for d in \
    b3/degroote \
    mesa/* \
    poly/* \
    hom; do

    cd $d

    ./gyre_ad < gyre_ad.in | tee gyre_ad.out

    cd -

done
        