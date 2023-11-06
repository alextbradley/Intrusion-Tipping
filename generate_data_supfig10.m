% Generate data for supplementary figure 10, showing how tidal velocities
% modulate the intrusion distance.

clear
addpath('functions');
tic

%% Parameters
timescale = 1; %day
u_inf = 0.01;  %m/s, far field velocity
u_tidal_dimensional = [0.001,0.005, 0.01,0.05, 0.1]; %m/s 

u_tidal   = u_tidal_dimensional/u_inf; % dimensionless tidal velocity

dT = [0.27, 0.3]; %dimensionless melt rate (=M in the manuscript)
S  = 0;           %dimensionless bedslope
F  = 0.25;        %Froude number
C  = 0.1;         %drag coefficient
                
%timestepping parameters
dt = 1e-3;
t_end = 100;
n_out = round(t_end/dt/20); %want this to be high enough that tidal oscillations are resolved
idx_out = round(linspace(0,t_end/dt,n_out));%unique(round(logspace(0,log10(sol(2).t_end /sol(2).dt), 100)));

%grid parameters
M = 100;
N = 100; 

for i = 1:2
    for j = 1:length(u_tidal)
        clear pp
        pp.u_tidal = u_tidal(j);
        pp.dT      = dT(i);  %=M in the manuscript
        pp.lambda  = 1/pp.dT;
        pp.F = F;
        pp.C = C;
        pp.S = S; 


        sol(i,j).pp = pp;
        sol(i,j).dt = dt;
        sol(i,j).t_end = t_end;
        sol(i,j).idx_out = idx_out;
        sol(i,j).t_plot = 2*t_end; %never plot

        %solve the system
        [sol(i,j).x,sol(i,j).t, sol(i,j).h, sol(i,j).h1, sol(i,j).u1, sol(i,j).Fr, sol(i,j).m, sol(i,j).dhdx,sol(i,j).dhdt,sol(i,j).idxmax] = evolve_shape_tidal(sol(i,j).pp.F, sol(i,j).pp.C, sol(i,j).pp.S, sol(i,j).pp.lambda,sol(i,j).pp.u_tidal,timescale, N, M, sol(i,j).dt, sol(i,j).t_end, sol(i,j).idx_out, sol(i,j).t_plot);


     end
end

% save the data
save('data-for-figures/supfigure10_data.mat', 'sol');
