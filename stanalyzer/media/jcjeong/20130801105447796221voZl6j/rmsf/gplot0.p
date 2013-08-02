set terminal png
set xlabel 'Atom Index'
set ylabel 'RMSF (Angstrom)'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130801105447796221voZl6j/rmsf/rmsf0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130801105447796221voZl6j/rmsf/rmsf_ouput_atom.dat" using 1:2 title "RMSF" with lines lw 3
