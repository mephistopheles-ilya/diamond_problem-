#!/bin/bash


./get --projections=400 --init_file=InitialModels/InitialModel5 --directory=init_examples1

./get --projections=400 --init_file=InitialModels/InitialModel5 --directory=init_examples2 --shift_poly3d=equations

./tran --directory_in=init_examples1 --directory_out=res_examples1 --projection=400 --shift=0 --points_in_contour=2000

./tran --directory_in=init_examples2 --directory_out=res_examples2 --projection=400 --shift=1 --points_in_contour=2000

./sd init_examples1 init_examples2 400 > sym_diff_compare.txt

./gc res_examples1 res_examples2 400 > gorizontal_compare.txt

