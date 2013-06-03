set terminal png
set xlabel 'Time (ps)'
set ylabel 'Angstrom'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130531165119115301bVYHWj/thickness_carbon/thickness_carbon0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130531165119115301bVYHWj/thickness_carbon/thickness_carbon.dat" using 1:2 title "phosphate thicknesss" with lines lw 3
