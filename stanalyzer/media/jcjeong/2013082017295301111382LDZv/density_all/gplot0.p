set terminal png
set encoding iso_8859_1
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/2013082017295301111382LDZv/density_all/density_all0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/2013082017295301111382LDZv/density_all/density_all.dat" using 1:2 title "Density" with lines lw 3
