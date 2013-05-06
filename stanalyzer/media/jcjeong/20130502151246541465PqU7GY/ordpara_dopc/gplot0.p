set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130502151246541465PqU7GY/ordpara_dopc/ordpara_dopc0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130502151246541465PqU7GY/ordpara_dopc/ordpara_depe_dcd710_2_pbs.dat" using 1:2 title "DOPC" with lines lw 3
