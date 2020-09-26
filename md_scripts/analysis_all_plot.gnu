# gnuplot file
# Jedd Bellamy-Carter, April 2020
# Requires files `rmsd.xvg`, `gyrate.xvg` and `area.xvg` in directory
# 
set terminal pngcairo enhanced font ",10"

set key left top
set grid

set output 'comp_rmsd.png'
set title 'RMSD of Production'
set xlabel 'Time (ps)'
set ylabel 'RMSD (nm)'
plot 'RT/rmsd.xvg' u 1:2 w l title 'RT','A/rmsd.xvg' u 1:2 w l title 'A','B/rmsd.xvg' u 1:2 w l title 'B','C/rmsd.xvg' u 1:2 w l title 'C'

set output 'comp_gyrate.png'
set title 'Gyration of Production'
set xlabel 'Time (ps)'
set ylabel 'R_g (nm)'
plot 'RT/gyrate.xvg' u 1:2 w l title 'RT','A/gyrate.xvg' u 1:2 w l title 'A','B/gyrate.xvg' u 1:2 w l title 'B','C/gyrate.xvg' u 1:2 w l title 'C'

set output 'comp_area.png'
set title 'Surface Area of Production'
set xlabel 'Time (ps)'
set ylabel 'Surface Area (nm^2)'
plot 'RT/area.xvg' u 1:2 w l title 'RT','A/area.xvg' u 1:2 w l title 'A','B/area.xvg' u 1:2 w l title 'B','C/area.xvg' u 1:2 w l title 'C'
