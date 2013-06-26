set data style linespoints
set border 31 lw 2.5
set title ''
set xlabel 'RMSD with respect to starting structure'
set ylabel 'Interaction Energy (kcal/mol)' 0.5,0.0
set ytics nomirror
set xtics nomirror
#set xrange [0:1]
#set yrange [0:1]
#set logscale x
#set logscale y

plot 'landscape.txt' using 2:4 title '' with linespoints 

pause -1 'Hit Return To Continue'
set style line	1 lt 1 lw 6 pt 6 ps 0.9  #red open circles
set style line	2 lt 3 lw 6 pt 7 ps 1.2  #blue closed circles
set style line	3 lt 2 lw 6 pt 4 ps 0.9  #green open squares
set style line	4 lt 4 lw 6 pt 5 ps 0.9  #purple closed squares
set style line	5 lt 5 lw 6 pt 8 ps 0.9  #cyan open triangles
set style line	6 lt 6 lw 6 pt 9 ps 0.9  #yellow closed triangles
set style line	7 lt 7 lw 6 pt 10 ps 0.9  #black open inverted triangles
set style line	8 lt 8 lw 6 pt 11 ps 0.9  #orange closed inverted triangles
set terminal postscript landscape color 22
set terminal postscript landscape solid
set terminal postscript landscape enhanced
set output 'plot.ps'

plot 'landscape.txt' using 2:4 title '' with linespoints ls 1 

