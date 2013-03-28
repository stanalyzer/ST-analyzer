set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130327133201811819Z6n5oX/density/density0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130327133201811819Z6n5oX/density/output.dat" using 1:2 title "Density" with lines lw 3
