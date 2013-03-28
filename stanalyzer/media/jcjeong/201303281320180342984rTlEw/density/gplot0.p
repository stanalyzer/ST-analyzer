set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201303281320180342984rTlEw/density/density0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201303281320180342984rTlEw/density/output.dat" using 1:2 title "Density" with lines lw 3
