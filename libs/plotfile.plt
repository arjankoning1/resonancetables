set title "Resonance Integral" font ",25"
set xlabel "C/E" font ",30"
set ylabel "#" font ",30" offset -2
set style line 1 linetype 1 linewidth 3 linecolor rgb "yellow"
set style line 2 linetype 1 linewidth 3 linecolor rgb "red"
set style line 3 linetype 1 linewidth 3 linecolor rgb "blue"
set style line 4 linetype 1 linewidth 3 linecolor rgb "green"
set style line 5 linetype 1 linewidth 3 linecolor rgb "black"
set term postscript enhanced font "Times-Roman,16" color
set output "bin.RI.eps"
plot \
"cendl3.2_bin.RI" u 1:2   title "CENDL-3.2" w histeps linestyle 1, \
"jeff3.3_bin.RI" u 1:2  title "JEFF3.3" w histeps linestyle 2, \
"endfb8.0_bin.RI" u 1:2   title "ENDFB8.0" w histeps linestyle 3, \
"jendl5.0_bin.RI" u 1:2   title "JENDL5.0" w histeps linestyle 4, \
"tendl.2023_bin.RI" u 1:2   title "TENDL-2023" w histeps linestyle 5, \
#   "dum" u 1:2 notitle
