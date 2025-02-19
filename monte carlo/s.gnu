set term png
set rmargin 5

set grid

set title "Convergence"
set title font "Arial, 12"

set xlabel 'N (Sample Size)'
set xtics font "Arial, 12"

set ylabel "σ"
set ytics font "Arial, 12"

set key right
set key font "Arial, 12"
set key spacing 2 R font "Arial,10"
set key box lt -1 lw 1

set output 'fig1.png'
plot 'res.dat' i 0 u 1:2 w l t "σ_{Real_1}",'res.dat' i 0 u 1:3 w l t "σ_{Estimated_1}" \
,'res.dat' i 0 u 1:4 w l t "σ_{Real_2}",'res.dat' i 0 u 1:5 w l t "σ_{Estimated_2}"

set output 'fig2.png'
plot 'res.dat' i 2 u 1:2:3 w yerrorbars t"Integral value"

set output '345.png'

set title 'σ / Approximated Multidimensional Integral'
plot 'res.dat' i 1 u 1:9 w l t"#1" ,'res.dat' i 1 u 1:10 w l t"#2",'res.dat' i 1 u 1:11 w l t"#3"

set output 'figMulti.png'
plot 'res.dat' i 2 u 1:2 w l t "Integral Value"