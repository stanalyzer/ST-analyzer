set terminal png
set multiplot layout 2, 1 title 'Area per lipid'
set tmargin 2
set title 'Top Membrane'
set xlabel 'Time (ps)'
set ylabel 'Area per Lipid (A^{2})'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130712152350339862RKA4JJ/lipid_per_area_voro/lipid_per_area_voro0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130712152350339862RKA4JJ/lipid_per_area_voro//home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130712152350339862RKA4JJ/lipid_per_area_voro/top_ave_voronoi_lipid_per_area.dat" using 1:2 title "['CHL1', 'POPC', 'POPI', 'POPI13', 'DOPC']" with lines lw 3, 
set title 'Bottom Membrane'
set xlabel 'Time (ps)'
set ylabel 'Area per Lipid (A^{2})'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130712152350339862RKA4JJ/lipid_per_area_voro/lipid_per_area_voro0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130712152350339862RKA4JJ/lipid_per_area_voro//home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130712152350339862RKA4JJ/lipid_per_area_voro/btm_ave_voronoi_lipid_per_area.dat" using 1:2 title "['CHL1', 'POPC', 'POPI', 'POPI13', 'DOPC']" with lines lw 3, 
