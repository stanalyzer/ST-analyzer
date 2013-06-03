set terminal png
set xlabel 'Time (ps)'
set ylabel 'RMSD (Angstrom)'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130603142759484993ivbYhC/rmsd/rmsd0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130603142759484993ivbYhC/rmsd/rmsd_ca_proa_ouput.dat" using 1:2 title "RMSD" with lines lw 3
