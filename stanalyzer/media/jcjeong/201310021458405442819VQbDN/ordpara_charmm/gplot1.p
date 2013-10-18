set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201310021458405442819VQbDN/ordpara_charmm/ordpara_charmm1.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201310021458405442819VQbDN/ordpara_charmm/ordpara_POPI.dat" using 1:2 title "POPI" with lines lw 3
