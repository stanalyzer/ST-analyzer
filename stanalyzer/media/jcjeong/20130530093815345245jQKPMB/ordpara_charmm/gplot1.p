set terminal png
set xlabel 'Carbon Index'
set ylabel 'S_CD'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130530093815345245jQKPMB/ordpara_charmm/ordpara_charmm1.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130530093815345245jQKPMB/ordpara_charmm/ordpara_DPPE_C3_300_561.dat" using 1:2 title "DOPC" with lines lw 3
