set data style linespoints
set border 31 lw 2.5
set title ''
set xlabel ''
set ylabel ''
set ytics nomirror
set xtics nomirror
#set xrange [0:1]
#set yrange [0:1]
#set logscale x
#set logscale y

plot 'gofr_caca_AA.dat' using 1:2 title 'columns 1:2' with linespoints, 'gofr_caca_AA.dat' using 1:3 title 'columns 1:3' with linespoints 

pause -1 'Hit Return To Continue'
set style line	1 lt 1 lw 6 pt 6 ps 0.9  #red open circles
set style line	2 lt 2 lw 6 pt 7 ps 1.2  #green closed circles
set style line	3 lt 3 lw 6 pt 4 ps 0.9  #blue open squares
set style line	4 lt 4 lw 6 pt 5 ps 0.9  #purple closed squares
set style line	5 lt 5 lw 6 pt 8 ps 0.9  #cyan open triangles
set style line	6 lt 6 lw 6 pt 9 ps 0.9  #yellow closed triangles
set style line	7 lt 7 lw 6 pt 10 ps 0.9  #black open inverted triangles
set style line	8 lt 8 lw 6 pt 11 ps 0.9  #orange closed inverted triangles
set terminal postscript landscape color 22
set terminal postscript landscape solid
set terminal postscript landscape enhanced
set output 'plot.ps'

plot 'gofr_caca_AA.dat' using 1:2 title 'columns 1:2' with linespoints ls 1, 'gofr_caca_AA.dat' using 1:3 title 'columns 1:3' with linespoints ls 2 

