set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130503145538097646hGrf7g/ordpara_dopc/ordpara_dopc0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130503145538097646hGrf7g/ordpara_dopc/ordpara_dppe_pbs_10.dat" using 1:2 title "DOPC" with lines lw 3
