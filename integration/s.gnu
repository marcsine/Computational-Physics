set term png

set output 'fig1.png'
set logscale x
set grid


plot 'p4-1920-res2.dat' u 1:2 w l t"res2 #1", 'p4-1920-res2.dat' u 1:3 w l t"res2 #2"

set yrange[40.9:42]

###################################

set output 'fig2.png'

plot 'p4-1920-res3.dat' u 1:2 w l t"res3 #1", 'p4-1920-res3.dat' u 1:3 w l t"res3 #2"


################################################

set output 'fig3.png'
set term pngcairo font ",10"
set format y "%.12f"
set multiplot layout 1, 2

set title 'Modified'
set yrange [41.668828320151:41.668828320156]
unset xtics

plot 'p4-1920-res4.dat' u 1:2 w l t'Trapezoidal',\
	'p4-1920-res4.dat' u 1:3 w l t'Simpson'

set yrange [40.9:41.8]
set title 'Unmodified'
unset format y

plot 'p4-1920-res3.dat' u 1:2 w l t'Trapezoidal', \
'p4-1920-res3.dat' u 1:3 w l t'Simpson'

unset multiplot