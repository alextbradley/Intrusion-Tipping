% Generate data for supplementary figure 8, showing the intrusion distance
% as a function of time for different slopes.

slopes = [-0.05, -0.01, 0, 0.01, 0.05];
cd = 1e-2; %drag coefficient

for i = 1:length(slopes)
    pp(i).F = 0.25;
    pp(i).dT = 0.30; %unbounded intrusion case in figure 2
    pp(i).lambda = 1/pp(i).dT;
    pp(i).S = slopes(i)/cd;
    pp(i).C    = 0.1;    %dimensionless interfacial drag

    sol(i).pp = pp(i);
    sol(i).dt = 1e-3;
    sol(i).t_end = 200;
    sol(i).t_plot = 1e5;
    sol(i).idx_out = unique(round(logspace(0,log10(sol(2).t_end /sol(2).dt), 300)));
    [sol(i).x,sol(i).t, sol(i).h, sol(i).h1, sol(i).u1, sol(i).Fr, sol(i).m, sol(i).dhdx,sol(i).dhdt,sol(i).idxmax] = evolve_shape(sol(i).pp.F, sol(i).pp.C, sol(i).pp.S, sol(i).pp.lambda, N, M, sol(i).dt, sol(i).t_end, sol(i).idx_out, sol(i).t_plot);

end

save('data-for-figures/supfigure8_data.mat', 'sol');