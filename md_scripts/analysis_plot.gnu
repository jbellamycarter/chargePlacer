# gnuplot file
# Jedd Bellamy-Carter, April 2020
# Requires files `rmsd.xvg`, `gyrate.xvg` and `area.xvg` in directory
# 
set terminal pngcairo font ",10"
set grid

set output 'rmsd.png'
set title 'RMSD of Production'
set xlabel 'Time (ps)'
set ylabel 'RMSD (nm)'
plot 'rmsd.xvg' using 1:2 with lines title ''

set output 'gyrate.png'
set title 'Gyration of Production'
set xlabel 'Time (ps)'
set ylabel 'R_g (nm)'
plot 'gyrate.xvg' using 1:2 with lines title ''

set output 'area.png'
set title 'Surface Area of Production'
set xlabel 'Time (ps)'
set ylabel 'Surface Area (nm^2)'
plot 'area.xvg' using 1:2 with lines title ''
