#FIGURA 1: Mètode de Gauss-Seidel

set term png
set output "fig1.png"

set grid # aixo posa una mena de graella, xarxa... no se com dir-ho

set xzeroaxis
set yzeroaxis

set title "Temperature convergence - Gauss-Seidel"
set xlabel "Iterations"
set ylabel "Temperature / ºC"

# set xrange [0:10000]
# set yrange [0:110]

plot "res.dat" index 0 u 1:2 w l t"Tint=6ºC" , \
""  index 1 u 1:2 with lines t"Tint=19ºC" , \
"" index 2 u 1:2 with lines t"Tint=320ºC"

#_______________________________________________________________________
#_______________________________________________________________________

#FIGURA 2: Mètode de Jacobi

set term png
set output "fig2.png"

set xzeroaxis
set yzeroaxis

set title "Temperature convergence - Jacobi"
set xlabel "Iterations"
set ylabel "Temperature / ºC"

#set xrange [0:10000]
#set yrange [10:50]

plot "res.dat" index 3 u 1:2 with lines t"Tint=6ºC" , \
"" index 4 u 1:2 with lines t"Tint=19ºC" , \
"" index 5 u 1:2 with lines t"Tint=320ºC"

#_______________________________________________________________________
#______________________________________________________________________

#FIGURA 3: Mètode de over-relaxation

set term png
set output "fig3.png"

set xzeroaxis
set yzeroaxis

set title "Temperature convergence - over-relaxation"
set xlabel "Iterations"
set ylabel "Temperature / ºC"

#set xrange [0:10000]
#set yrange [10:50]

plot "res.dat" index 6 u 1:2 with lines t"Tint=6ºC" , \
""  index 7 u 1:2 with lines t"Tint=19ºC" , \
"" index 8 u 1:2 with lines t"Tint=320ºC"

#_______________________________________________________________________
#______________________________________________________________________


#FIGURA 4:Distribució de temperatures

set term png
set output "fig4.png"

set xzeroaxis
set yzeroaxis
set zzeroaxis


set title "Final Temperature distribution"
set xlabel "x / cm"
set ylabel "y / cm"

# set xrange [0:40]
# set yrange [0:27.5]
# set zrange [10:35]

set pm3d
set view 0.0

splot "res2.dat" u 1:2:3 w l t"Temperature / ºC

