set terminal png
set xlabel 'Time (ps)'
set ylabel 'Angstrom'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201308011636126283543pgoOr/thickness_custom/thickness_custom0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201308011636126283543pgoOr/thickness_custom/thickness_custom_0.dat" using 1:2 title "custom thicknesss" with lines lw 3
