#!/bin/bash

./stress_test.py -c 10 -n ad_matrix_solver -v BANDED,CYCLIC,ROWPP gyre.in.template_ad stress.ad_matrix_solver.csv

./stress_test.py -c 10 -t 1,4,8 -n ad_matrix_solver -v BANDED,CYCLIC,ROWPP gyre.in.template_ad stress.ad_matrix_solver.multi.csv
./stress_test.py -c 10 -t 1,4,8 -n nad_matrix_solver -v BANDED,CYCLIC,ROWPP gyre.in.template_nad stress.nad_matrix_solver.multi.csv
