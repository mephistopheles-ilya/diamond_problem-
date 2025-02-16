#!/bin/bash

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=0 --parametr_of_shift=0.00001 --distribution=normal --sigma=0.00001
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model0.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.00001 --distribution=normal --sigma=0.00001
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model1.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.00005 --distribution=normal --sigma=0.00005
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model2.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.0001 --distribution=normal --sigma=0.0001
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model3.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.0005 --distribution=normal --sigma=0.0005
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model4.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.001 --distribution=normal --sigma=0.001
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model5.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.005 --distribution=normal --sigma=0.005
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model6.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.01 --distribution=normal --sigma=0.01
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model7.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.05 --distribution=normal --sigma=0.05
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model8.pdf

./a.out --initial_model=InitialModels/InitialModel7 --is_change_shift=1 --parametr_of_shift=0.1 --distribution=normal --sigma=0.1
gnuplot -e "filename_out='photos/photo.png'; filename_in='gnu.txt'" gnu_plot.gp
convert photos/photo.png photos/model9.pdf

