set logscale y
set format y "%e";
set title "Residuals"
set ylabel 'Residual'
set xlabel  'Iteration'
f(x) = 1e-05
g(x) = 1e-06
plot"< cat ../log.oildensityFoam | grep 'Solving for p'    | cut -d' ' -f9 | tr -d ','" title 'p' with lines, \
"< cat ../log.oildensityFoam | grep 'Solving for psi'    | cut -d' ' -f9 | tr -d ','" title 'psi' with lines, \
"< cat ../log.oildensityFoam | grep 'Solving for lamda'    | cut -d' ' -f9 | tr -d ','" title 'lamda' with lines, \
"< cat ../log.oildensityFoam | grep 'Solving for T'    | cut -d' ' -f9 | tr -d ','" title 'T' with lines, \
"< cat ../log.oildensityFoam | grep 'Solving for epsilon'    | cut -d' ' -f9 | tr -d ','" title 'epsilon' with lines, \
"< cat ../log.oildensityFoam | grep 'Solving for k'    | cut -d' ' -f9 | tr -d ','" title 'k' with lines, \
f(x) title '1e-05' with lines, \
g(x) title '1e-06' with lines
pause 1
reread

/*
"< cat ../log.oildensityFoam | grep 'Solving for Ux'    | cut -d' ' -f9 | tr -d ','" title 'Ux' with lines, \
"< cat ../log.oildensityFoam | grep 'Solving for Uy'    | cut -d' ' -f9 | tr -d ','" title 'Uy' with lines, \
"< cat ../log.oildensityFoam | grep 'Solving for c'    | cut -d' ' -f9 | tr -d ','" title 'c' with lines, \

*/
