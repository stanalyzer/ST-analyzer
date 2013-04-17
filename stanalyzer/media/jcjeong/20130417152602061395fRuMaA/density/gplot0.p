set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417152602061395fRuMaA/density/density0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417152602061395fRuMaA/density/density_custom1.dat" using 1:2 title "Density" with lines lw 3
