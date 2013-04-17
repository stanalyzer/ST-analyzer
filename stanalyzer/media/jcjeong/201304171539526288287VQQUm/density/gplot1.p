set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201304171539526288287VQQUm/density/density1.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201304171539526288287VQQUm/density/density_custom1.dat" using 1:2 title "Density" with lines lw 3
