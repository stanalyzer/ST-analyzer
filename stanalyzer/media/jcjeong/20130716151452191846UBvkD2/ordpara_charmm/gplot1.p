set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716151452191846UBvkD2/ordpara_charmm/ordpara_charmm1.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716151452191846UBvkD2/ordpara_charmm/ordpara_POPI_1.dat" using 1:2 title "DOPC" with lines lw 3
