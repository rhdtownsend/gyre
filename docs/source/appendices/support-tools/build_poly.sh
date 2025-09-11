#!/bin/sh

$GYRE_DIR/bin/build_poly --n_poly=3 --dz=0.01 poly.simple.h5
$GYRE_DIR/bin/build_poly --n_poly=3,1.5 --z_b=1.4 --Delta_b=-0.5 --dz=0.01 poly.composite.h5
