set term wxt 1
set xlabel 'Time'
set ylabel 'p: x=Lx/2'
#set yrange [0.009:0.011]
#plot 'postProcessing/probes/0/p' u 1:2
#rep 'postProcessing/probes/0/p' u 1:3
plot 'postProcessing/probes/0/p' u 1:4
rep 'postProcessing/probes/0/p' u 1:5
#rep 'postProcessing/probes/0/p' u 1:6
#rep 'postProcessing/probes/0/p' u 1:7

set term wxt 2
set xlabel 'Time'
set ylabel 'p: x=Lx'
#set yrange [0.009:0.011]
#plot 'postProcessing/probes/0/p' u 1:2
#rep 'postProcessing/probes/0/p' u 1:3
#plot 'postProcessing/probes/0/p' u 1:4
#rep 'postProcessing/probes/0/p' u 1:5
plot 'postProcessing/probes/0/p' u 1:6
rep 'postProcessing/probes/0/p' u 1:7

set term wxt 3
set xlabel 'Time'
set ylabel 'Ux near wall'
#set yrange [0.85:0.86]
#plot 'postProcessing/probes/0/tmp' u 1:2
#rep 'postProcessing/probes/0/tmp' u 1:5
plot 'postProcessing/probes/0/tmp' u 1:8
#rep 'postProcessing/probes/0/tmp' u 1:11
rep 'postProcessing/probes/0/tmp' u 1:14
#rep 'postProcessing/probes/0/tmp' u 1:17

set term wxt 4
set xlabel 'Time'
set ylabel 'Ux middle'
#set yrange [0.85:0.86]
#plot 'postProcessing/probes/0/tmp' u 1:2
#plot 'postProcessing/probes/0/tmp' u 1:5
#rep 'postProcessing/probes/0/tmp' u 1:8
plot 'postProcessing/probes/0/tmp' u 1:11
#rep 'postProcessing/probes/0/tmp' u 1:14
rep 'postProcessing/probes/0/tmp' u 1:17

set term wxt 5
set xlabel 'Time'
set ylabel 'Uy near wall'
#set yrange [-1.e-5:1.e-5]
#plot 'postProcessing/probes/0/tmp' u 1:3
#rep 'postProcessing/probes/0/tmp' u 1:6
plot 'postProcessing/probes/0/tmp' u 1:9
#rep 'postProcessing/probes/0/tmp' u 1:12
rep 'postProcessing/probes/0/tmp' u 1:15
#rep 'postProcessing/probes/0/tmp' u 1:18

set term wxt 6
set xlabel 'Time'
set ylabel 'Uy middle'
#set yrange [-1.e-5:1.e-5]
#plot 'postProcessing/probes/0/tmp' u 1:3
#plot 'postProcessing/probes/0/tmp' u 1:6
#rep 'postProcessing/probes/0/tmp' u 1:9
plot 'postProcessing/probes/0/tmp' u 1:12
#rep 'postProcessing/probes/0/tmp' u 1:15
rep 'postProcessing/probes/0/tmp' u 1:18

pause -1