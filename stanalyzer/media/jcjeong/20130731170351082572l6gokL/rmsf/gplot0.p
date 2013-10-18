set terminal png
set xlabel 'Atom Index'
set ylabel 'RMSF (Angstrom)'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130731170351082572l6gokL/rmsf/rmsf0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130731170351082572l6gokL/rmsf/rmsf_ouput.dat" using 1:2 title "RMSF" with lines lw 3
