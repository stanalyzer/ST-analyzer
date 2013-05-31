set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130531105901714004RIg0oP/ordpara_charmm/ordpara_charmm2.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130531105901714004RIg0oP/ordpara_charmm/ordpara_DPPE_C2_150_750.dat" using 1:2 title "DOPC" with lines lw 3
