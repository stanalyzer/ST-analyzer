set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201305031402587358458S7iHd/ordpara_dopc/ordpara_dopc0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201305031402587358458S7iHd/ordpara_dopc/ordpara_dopc.dat" using 1:2 title "DOPC" with lines lw 3
