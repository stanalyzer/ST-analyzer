set terminal png
set xlabel 'Time (ps)'
set ylabel 'Angstrom'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201307161439244602537QP8xG/thickness_carbon/thickness_carbon0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201307161439244602537QP8xG/thickness_carbon/thickness_carbon.dat" using 1:2 title "carbon_based thicknesss" with lines lw 3
