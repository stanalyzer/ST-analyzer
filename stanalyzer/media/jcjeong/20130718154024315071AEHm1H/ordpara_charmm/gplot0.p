set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130718154024315071AEHm1H/ordpara_charmm/ordpara_charmm0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130718154024315071AEHm1H/ordpara_charmm/ordpara_dppe_0.dat" using 1:2 title "DPPE" with lines lw 3
