set terminal pngcairo size 800,600 enhanced font 'Arial,10'

set output filename_out

set view 0, 0  # Azimuth = 0, Elevation = 0

set xlabel 'X-axis'
set ylabel 'Y-axis'
set zlabel 'Z-axis'

set grid
set box

set title '3D Model View from Positive Z-axis'

splot filename_in with lines lw 2 title '3D Model'

unset output
