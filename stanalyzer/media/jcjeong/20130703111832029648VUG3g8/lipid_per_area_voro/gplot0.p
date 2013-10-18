set terminal png
set xlabel 'Time (ps)'
set ylabel 'Area per Lipid (A^2)'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130703111832029648VUG3g8/lipid_per_area_voro/lipid_per_area_voro0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130703111832029648VUG3g8/lipid_per_area_voro/voronoi_lipid_per_area.dat" using 1:2 title "custom thicknesss" with lines lw 3
