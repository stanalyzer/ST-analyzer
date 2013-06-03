set terminal png
set xlabel 'Time (ps)'
set ylabel 'RMSD (Angstrom)'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201306031414594993654l7hwq/rmsd/rmsd0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201306031414594993654l7hwq/rmsd/rmsd_ouput.dat" using 1:2 title "RMSD" with lines lw 3
