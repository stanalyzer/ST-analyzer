set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130503152417289013nLNC4C/ordpara_dopc/ordpara_dopc0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130503152417289013nLNC4C/ordpara_dopc/ordpara_dxpe_151_pbs_20.dat" using 1:2 title "DOPC" with lines lw 3
