set term png
set output 'fig1.png'

set grid

set title 'Big oscillations'
set title font "Arial, 12"

set xlabel 't / s'
set xtics font "Arial,11"

set ylabel 'dθ/dt / (rad/s)'
set ytics font "Arial,11"

set key out horiz top right
set key font "Arial, 9"
set key spacing 1.5 R font "Arial,10"
set key box lt -1 lw 1
set key spacing 1 R

plot 'res.dat' i 0 u 1:3 w l t 'Euler','res.dat' i 1 u 1:3 w l t 'Improved Euler'


##############

set output 'fig2.png'

set title 'Big oscillations, phase space'
set title font "Arial, 12"

set xlabel 'θ / rad'
set xtics font "Arial,11"

set ylabel 'dθ/dt / (rad/s)'
set ytics font "Arial,11"


plot 'res.dat' i 0 u 2:3 w l t 'Euler','res.dat' i 1 u 2:3 w l t 'Improved Euler'

########################

set output 'fig3.png'

set title 'Small oscillations, phase space'
set title font "Arial, 12"

set xlabel 't / s'
set xtics font "Arial,11"

set ylabel 'dθ/dt / (rad/s)'
set ytics font "Arial,11"

plot 'res.dat' i 2 u 1:3 w l t 'Euler','res.dat' i 3 u 1:3 w l t 'Improved Euler'

########################

set output 'fig4.png'

set title 'Potential and Total Energy vs Time'
set title font "Arial, 12"

set xlabel 't / s'
set xtics font "Arial,11"

set ylabel 'Energies / J'
set ytics font "Arial,11"

plot 'res.dat' i 4 u 1:5 w l t 'U Euler θ(0)=1',\
'res.dat' i 4 u 1:6 w l t 'E Euler θ(0)=1',\
'res.dat' i 5 u 1:5 w l t 'U Euler Imp. θ(0)=1',\
'res.dat' i 5 u 1:6 w l t 'E Euler Imp. θ(0)=1',\
'res.dat' i 6 u 1:5 w l t 'U Euler θ(0)=π-0.02',\
'res.dat' i 6 u 1:6 w l t 'E Euler θ(0)=π-0.02',\
'res.dat' i 7 u 1:5 w l t 'U Euler Imp. θ(0)=π-0.02',\
'res.dat' i 7 u 1:6 w l t 'E Euler Imp. θ(0)=π-0.02'

######################

set output 'fig5.png

set title 'Phase space'
set title font "Arial, 12"

set xlabel 'θ / rad'
set xtics font "Arial,11"

set ylabel 'dθ/dt / (rad/s)'
set ytics font "Arial,11"


plot 'res.dat' i 8 u 2:3 w d t '+', 'res.dat' i 9 u 2:3 w d t '-'

#######################

set output 'fig6.png' 

set title 'Convergence'
set title font "Arial, 12"

set xlabel 't / s'
set xtics font "Arial,11"

set ylabel 'Total Energy'
set ytics font "Arial,11"

plot 'res.dat' i 10 u 1:6 w l t '600 time steps',\
'res.dat' i 11 u 1:6 w l t '1300 time steps',\
'res.dat' i 12 u 1:6 w l t '2600 time steps',\
'res.dat' i 13 u 1:6 w l t '15000 time steps'


