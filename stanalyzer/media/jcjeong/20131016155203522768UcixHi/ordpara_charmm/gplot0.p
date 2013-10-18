set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20131016155203522768UcixHi/ordpara_charmm/ordpara_charmm0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20131016155203522768UcixHi/ordpara_charmm/Ecoli_pe27_1_ordpara_charmm_DXPE_C2.dat" using 1:2 title "DXPE" with lines lw 3
