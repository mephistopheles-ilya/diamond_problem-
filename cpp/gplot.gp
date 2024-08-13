#set terminal wxt size 800,800 enhanced font 'Arial,12' persist
set title "Plot of 2D Points"
set xlabel "X Axis"
set ylabel "Y Axis"
set grid
set style line 1 lc 'red' lw 0.7 pt 7 ps 0.5
set style line 2 lc 'black' lw 0.7 pt 7 ps 0.5 
set style line 3 lc 'blue' dt(3, 8)  lw 2 pt 7 ps 0.1
set key left top
#plot [-0.5:1.25] [-0.5:1.25] 'in_triangle.txt' with linespoints ls 1 title "File 1", 'out_triangle.txt' with linespoints ls 2 title "File 2", 'triangle_intersection.txt' with polygons lc 'blue' fs solid 0.5 title "File 3"
# plot 'in_pentagon.txt' with linespoints ls 1 title "File 1", 'out_pentagon.txt' with linespoints ls 2 title "File 2"
plot [-1.75:1.75] [-1.75:1.75] 'in_decagon.txt' with linespoints ls 1 title "File 1", 'out_decagon.txt' with linespoints ls 2 title "File 2" , 'decagon_intersection.txt' with polygons lc 'blue' fs solid 0.5 title "File 3"
pause -1 "Hit any key to continue"
