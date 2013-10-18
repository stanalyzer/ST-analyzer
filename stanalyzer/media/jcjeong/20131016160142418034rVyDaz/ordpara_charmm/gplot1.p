set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20131016160142418034rVyDaz/ordpara_charmm/ordpara_charmm1.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20131016160142418034rVyDaz/ordpara_charmm/ordpara_charmm_DXPE_C3.dat" using 1:2 title "DXPE" with lines lw 3
