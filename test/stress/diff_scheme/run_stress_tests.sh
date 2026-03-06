#!/bin/bash

./stress_test.py -c 10 -n diff_scheme -v COLLOC_GL2,COLLOC_GL4,COLLOC_GL6,MAGNUS_GL2,MAGNUS_GL4,MAGNUS_GL6,TRAPZ,MIRK \
		 gyre.in.template stress.ad_matrix_solver.csv
