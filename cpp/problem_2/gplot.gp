set title "Plot of 2D Points"
set xlabel "X Axis"
set ylabel "Y Axis"
set grid
set style line 1 lc 'red' lw 0.7 pt 7 ps 0.7
set style line 2 lc 'green' lw 0.7 pt 7 ps 0.7 
set key left top

plot [-5:5] [-7:-1] 'points2d_before.txt' with linespoints ls 1 title "with wrong part", 'points2d_after.txt' with linespoints ls 2 title "wrong part deleted" 

pause -1 "Hit any key to constineu"
