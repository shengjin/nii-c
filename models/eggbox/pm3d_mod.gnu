reset

#load 'brbg.pal'
#load 'rdbu.pal'
#load 'rdgy.pal'
#load 'pubu.pal'
#load 'blues.pal'
#load 'bugn.pal'
#load 'greens.pal'
load 'bentcoolwarm.pal'
#load 'moreland.pal'  




set term pdf size 2.4,2 linewidth 1.5 font "NotoSerif,10"
set output "output.mod.pdf"

#set border 4095

set xrange [0:50]
set yrange [0:50]

set xtics out 0,10,60  offset 0,0.5
set ytics out 0,10,50  offset 0.5,0
set mxtics 10
set mytics 10

#set palette @viridis
#set palette @plasma
#set palette @magma
#set palette @MATLAB
#set palette @inferno
#set palette @rdbu 
#set palette 

set colorbox 

set pm3d map
set cbtics 0,1,2 font "NotoSerif, 9" offset -0.8,0
set cbrange [0:2]
set mcbtics 10
#set style data pm3d
#set pm3d depthorder

#set view projection xy

set xlabel "{/:Italic X}"  font "NotoSerif, 14" offset 0,0.8 
set ylabel "{/:Italic Y}"  font "NotoSerif, 14" offset 1.0,0 
set cblabel "{/:Italic Posterior}"  font "NotoSerif, 10" offset -0.6,-0.20 

plot "xy_fit_mod.dat" u 1:2:4 w p linecolor palette pt 5 ps 0.1 title ""


