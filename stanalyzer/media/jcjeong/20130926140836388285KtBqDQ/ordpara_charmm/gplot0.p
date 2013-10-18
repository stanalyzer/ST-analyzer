set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130926140836388285KtBqDQ/ordpara_charmm/ordpara_charmm0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130926140836388285KtBqDQ/ordpara_charmm/ordpara_charmm_0.dat" using 1:2 title "POPC" with lines lw 3
