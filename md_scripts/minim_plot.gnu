# gnuplot file
# Jedd Bellamy-Carter, April 2020
# Requires file "min_energy.xvg" in directory
# 
set terminal pngcairo enhanced font ",10"
set grid
set output 'min_energy.png'
set title 'Energy Minimisation'
set xlabel 'Steps (n)'
set ylabel 'Potential Energy (kJ/mol)'
plot 'min_energy.xvg' using 1:2 with lines title ''
