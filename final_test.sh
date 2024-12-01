#!/bin/bash

echo "Get projections"

./get --projections=400 --init_file=InitialModels/InitialModel5 --directory=init_examples1

./get --projections=400 --init_file=InitialModels/InitialModel5 --directory=init_examples2 --shift_poly3d=equations --distribution=uniform --is_change_slope=1 --parametr_of_slope=0.0001 --is_change_azimuth=1 --parametr_of_azimuth=0.0001



./a.out --initial_model=InitialModels/InitialModel5 --output=tmp.ply
./a.out --initial_model=InitialModels/InitialModel5 --output=tmp_changed.ply --distribution=uniform --is_change_slope=1 --parametr_of_slope=0.0001  --is_change_azimuth=1 --parametr_of_azimuth=0.0001


#echo "Spoil projections"

#time ./tran --directory_in=init_examples1 --directory_out=res_examples1 --projection=400 --shift=0 --points_in_contour=2000

#time ./tran --directory_in=init_examples2 --directory_out=res_examples2 --projection=400 --shift=0 --points_in_contour=2000 --sz_shift=0.0001

echo "Compare prijections"

./sd init_examples1 init_examples2 400 > sym_diff_compare.txt
#time ./sd res_examples1 res_examples2 400 > sym_diff_compare_2.txt

#time ./gc res_examples1 res_examples2 400 > gorizontal_compare.txt

