#set terminal wxt size 800,800 enhanced font 'Arial,12' persist
set title "Plot of 2D Points"
set xlabel "X Axis"
set ylabel "Y Axis"
set grid
set style line 1 lc 'red' lw 0.7 pt 7 ps 0.5
set style line 2 lc 'green' lw 0.7 pt 7 ps 0.5 
#set style line 3 lc 'blue' dt(3, 8)  lw 2 pt 7 ps 0.1
set style line 3 lc 'blue' lw 0.7 pt 7 ps 0.5
set key left top

plot [-5:5] [-7:-1] 'poly2d-cont000' with linespoints ls 1 title "input figure", 'out_poly2d-cont000' with linespoints ls 2 title "output figure" #, 'poly_intersection' with linespoints ls 3 title "dif"


#plot [-5:5] [-7:-1] 'poly2d-cont000' with linespoints ls 1 title "input figure", 'poly_intersection' with polygons lc 'blue' fs solid 0.5 title "difference"

pause -1 "Hit any key to constineu"
