set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130821110540708466OSd3hG/ordpara_charmm/ordpara_charmm0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130821110540708466OSd3hG/ordpara_charmm/ordpara_charmm_0.dat" using 1:2 title "DOPC" with lines lw 3
