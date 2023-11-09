% Generate data for supplementary figure 8, showing the intrusion distance
% as a function of time for different slopes.
addpath('functions');
Sc = 0.1678; %for no melting with F = 0.25, C = 0.1 (i.e. for slopes beyond this, we have unbounded intrusion at the first timestep)
cd = 1e-2; %drag coefficient

slopes_unscaled = [-0.05, -0.01,-0.005,-0.001, 0, 0.0008, 0.0016]; 
slopes_scaled = slopes_unscaled/cd;

intrusion_lengths = cell(1,length(slopes_scaled));
times   =cell(1,length(slopes_scaled));

for i = 1:length(slopes_scaled)
    N = 200; M = 50; %number of blocks and grid size
    pp(i).F = 0.25;
    pp(i).dT = 0.27; %bounded intrusion case in figure 2
    pp(i).lambda = 1/pp(i).dT;
    pp(i).S = slopes_scaled(i);
    pp(i).C    = 0.1;    %dimensionless interfacial drag

    sol(i).pp = pp(i);
    sol(i).dt = 1e-3;
    sol(i).t_end = 100;
    sol(i).t_plot = 1e5;
    sol(i).idx_out = unique(round(logspace(0,log10(sol(i).t_end /sol(i).dt), 100)));
    [sol(i).x,sol(i).t, sol(i).h, sol(i).h1, sol(i).u1, sol(i).Fr, sol(i).m, sol(i).dhdx,sol(i).dhdt,sol(i).idxmax] = evolve_shape(sol(i).pp.F, sol(i).pp.C, sol(i).pp.S, sol(i).pp.lambda, N, M, sol(i).dt, sol(i).t_end, sol(i).idx_out, sol(i).t_plot);

      %get the intrusion length
        intrusion_length = nan(1,length(sol(i).t));
        t = nan(1,length(sol(i).t));
        for iT = 1:(length(sol(i).t)-1)
            idx = find(sol(i).h(iT,:)-sol(i).h1(iT,:) == 0, 1, 'First');
            intrusion_length(iT) = -sol(i).x(idx);
            t(iT) = sol(i).t(iT);
        end

        intrusion_lengths{i} = intrusion_length;
        times{i} = t;

        %save the intrusion length
        save('data-for-figures/supfigure8_data_.mat', 'intrusion_lengths', 'times');
end


