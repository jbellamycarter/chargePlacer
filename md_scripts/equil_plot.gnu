# gnuplot file
# Jedd Bellamy-Carter, April 2020
# Requires file `eq_energy.xvg` in directory
# 
set terminal pngcairo font ",10"
set grid

set output 'eq_temp.png'
set title 'Temperature of Equilibration'
set xlabel 'Time (ps)'
set ylabel 'Temperature (K)'
plot 'eq_energy.xvg' using 1:4 with lines title ''

set output 'eq_energy.png'
set title 'Energy of Equilibration'
set xlabel 'Time (ps)'
set ylabel 'Energy (kJ/mol)'
plot 'eq_energy.xvg' using 1:2 with lines title 'Coulomb', '' using 1:3 with lines title 'Potential'
