set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417154434003921YQNbiV/density/density3.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417154434003921YQNbiV/density/density_water.dat" using 1:2 title "Density" with lines lw 3
