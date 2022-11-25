% Generate the data used in figure 2 of the manuscript.

% Generate outputs for the evolution of the channel shape. Outputs are for
% for two different values of dT. The lower value has a bounded intrusion,
% while the higher values (more melting) has an unbounded intrusion.  Here
% we have only a flat bed (S = 0).
%
% AB 21/11/22 MIT License
%% Preliminaries
clear
addpath('functions');
tic
%% Parameters 
% set one
pp(1).F    = 0.25;   %froude number on entry Fr_infty = Q_{infty} / sqrt(g' * H_{infty}^3)
pp(1).C    = 0.1;    %dimensionless interfacial drag
pp(1).S    = 0;	 	 %bed slope
pp(1).dT   = 0.27;   %dimensionless ocean temperature
pp(1).lambda = 1/pp(1).dT;

%set two; 
pp(2)      = pp(1);  %make a copy
pp(2).dT   = 0.30; %adjust dT
pp(2).lambda = 1/pp(2).dT;

% seed the solutions
sol(1).pp = pp(1);
sol(2).pp = pp(2);
%% Grid and timestepping
N  	 = 100;     %number of grid points inside each block
M 	 = 20; 		%number of blocks

%timestepping
sol(1).dt 	 = 1e-3; %solution timestep
sol(1).t_end   = 100;    %end time for solution 1
sol(1).t_plot  = 1000;   %when to plot solution (make larger than t_end to suppress)
sol(1).t_store = 0.1;  %when to store solution
sol(1).idx_out = unique(round(logspace(0,log10(sol(1).t_end /sol(1).dt), 100)));

sol(2).dt 	 = 1e-3; %solution timestep
sol(2).t_end   = 100;    %end time for solution 2
sol(2).t_plot  = 1000; %when to plot solution 
%sol(2).t_store = 0.1;  %when to store solution
sol(2).idx_out = unique(round(logspace(0,log10(sol(2).t_end /sol(2).dt), 100)));
%% Get the solutions
[sol(1).x,sol(1).t, sol(1).h, sol(1).h1, sol(1).u1, sol(1).Fr, sol(1).m, sol(1).dhdx,sol(1).dhdt,sol(1).idxmax] = evolve_shape(sol(1).pp.F, sol(1).pp.C, sol(1).pp.S, sol(1).pp.lambda, N, M, sol(1).dt, sol(1).t_end, sol(1).idx_out, sol(1).t_plot);
[sol(2).x,sol(2).t, sol(2).h, sol(2).h1, sol(2).u1, sol(2).Fr, sol(2).m, sol(2).dhdx,sol(2).dhdt,sol(2).idxmax] = evolve_shape(sol(2).pp.F, sol(2).pp.C, sol(2).pp.S, sol(2).pp.lambda, N, M, sol(2).dt, sol(2).t_end, sol(2).idx_out, sol(2).t_plot);

%% Save the solution
save('data-for-figures/figure1_data.mat', 'sol');
toc
