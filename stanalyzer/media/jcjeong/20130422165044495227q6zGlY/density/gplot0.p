set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130422165044495227q6zGlY/density/density0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130422165044495227q6zGlY/density/density_water.dat" using 1:2 title "Density" with lines lw 3
