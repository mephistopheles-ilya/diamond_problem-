#!/bin/bash

#echo "Get projections"

#./get.out --projections=400 --init_file=InitialModels/InitialModel5 --directory=init_examples1

#./get.out --projections=400 --init_file=InitialModels/InitialModel5 --directory=init_examples2 --shift_poly3d=equations --distribution=uniform --is_change_slope=1 --parametr_of_slope=0.0001 --is_change_azimuth=1 --parametr_of_azimuth=0.0001



#./splot.out --initial_model=InitialModels/InitialModel5 --output=tmp.ply
#./splot.out --initial_model=InitialModels/InitialModel5 --output=tmp_changed.ply --distribution=uniform --is_change_slope=1 --parametr_of_slope=0.0001  --is_change_azimuth=1 --parametr_of_azimuth=0.0001


#echo "Spoil projections"

#time ./tran.out --directory_in=init_examples1 --directory_out=res_examples1 --projection=400 --shift=0 --points_in_contour=2000

#time ./tran.out --directory_in=init_examples2 --directory_out=res_examples2 --projection=400 --shift=0 --points_in_contour=2000 --sz_shift=0.0001

#echo "Compare prijections"

#./sym_diff.out init_examples1 init_examples2 400 > sym_diff_compare.txt
#time ./sym_diff.out res_examples1 res_examples2 400 > sym_diff_compare_2.txt

#time ./gorizontal_compare.out res_examples1 res_examples2 400 > gorizontal_compare.txt


./get.out --projections=400 --init_file=InitialModels/InitialModel5 --directory=init_examples
./tran.out --projections=400 --gr_method=ceils --dist_points=0.008 --shift=0
./tran.out --projections=400 --gr_method=ceils --dist_points=0.008 --shift=1 --sz_shift=0.005 --directory_out=res_examples1
./tran.out --projections=400 --gr_method=ceils --dist_points=0.008 --shift=1 --sz_shift=0.01 --directory_out=res_examples2

./txt_to_obj.out --file_txt=InitialModels/InitialModel5 --file_obj=model_00001.obj
./txt_to_con.out --name_of_con_file=model_00001_grid.con1 --directory=/home/ilya/diamond_problem-/res_examples --type_of_file=con1 --model_in_txt=InitialModels/InitialModel5
./txt_to_con.out --name_of_con_file=model_00001_grid.con2 --directory=/home/ilya/diamond_problem-/res_examples --type_of_file=con2 --model_in_txt=InitialModels/InitialModel5
./txt_to_con.out --name_of_con_file=model_00001_grid.con3 --directory=/home/ilya/diamond_problem-/res_examples --type_of_file=con3 --model_in_txt=InitialModels/InitialModel5
./txt_to_con.out --name_of_con_file=model_00001_grid_shift_0005.con1 --directory=/home/ilya/diamond_problem-/res_examples1 --type_of_file=con1 --model_in_txt=InitialModels/InitialModel5
./txt_to_con.out --name_of_con_file=model_00001_grid_shift_0005.con2 --directory=/home/ilya/diamond_problem-/res_examples1 --type_of_file=con2 --model_in_txt=InitialModels/InitialModel5
./txt_to_con.out --name_of_con_file=model_00001_grid_shift_0005.con3 --directory=/home/ilya/diamond_problem-/res_examples1 --type_of_file=con3 --model_in_txt=InitialModels/InitialModel5
./txt_to_con.out --name_of_con_file=model_00001_grid_shift_001.con1 --directory=/home/ilya/diamond_problem-/res_examples2 --type_of_file=con1 --model_in_txt=InitialModels/InitialModel5
./txt_to_con.out --name_of_con_file=model_00001_grid_shift_001.con2 --directory=/home/ilya/diamond_problem-/res_examples2 --type_of_file=con2 --model_in_txt=InitialModels/InitialModel5
./txt_to_con.out --name_of_con_file=model_00001_grid_shift_001.con3 --directory=/home/ilya/diamond_problem-/res_examples2 --type_of_file=con3 --model_in_txt=InitialModels/InitialModel5

