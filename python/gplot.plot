set terminal wxt size 800,800 enhanced font 'Arial,12' persist
set title "Plot of 2D Points"
set xlabel "X Axis"
set ylabel "Y Axis"
set grid
set style line 1 lc 'red' pt 7 ps 0.5
set style line 2 lc 'blue' pt 5 ps 0.5
set key left top
# plot 'in_triangle.txt' with linespoints ls 1 title "File 1", 'out_triangle.txt' with linespoints ls 2 title "File 2"
# plot 'in_pentagon.txt' with linespoints ls 1 title "File 1", 'out_pentagon.txt' with linespoints ls 2 title "File 2"
plot [-1.75:1.75] [-1.75:1.75] 'in_decagon.txt' with linespoints ls 1 title "File 1", 'out_decagon.txt' with linespoints ls 2 title "File 2"
