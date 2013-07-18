set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201307181517130444484eOOTa/ordpara_charmm/ordpara_charmm0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201307181517130444484eOOTa/ordpara_charmm/ordpara_popc_0.dat" using 1:2 title "DOPC" with lines lw 3
