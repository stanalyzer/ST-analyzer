set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130503145934485007gT0IW0/ordpara_dopc/ordpara_dopc0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130503145934485007gT0IW0/ordpara_dopc/ordpara_dppe_400.dat" using 1:2 title "DOPC" with lines lw 3
