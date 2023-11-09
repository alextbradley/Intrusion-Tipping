function [x,t_out, h_out, h1_out, u1_out, Fr_out, melt_out, dhdx_out, dhdt_out,idxmax] = evolve_shape_tidal(Fr_infty, C, S, lambda,u_tidal,timescale, N, M, dt, t_end, idx_out, t_plot)
% Evolve the shape of the channel. 
%
%% Inputs
% Fr_infty  : 	(scalar) froude number on entry Fr_infty = Q_{infty} / sqrt(g' * H_{infty}^3)
% C         :   (scalar) reduced interfacial drag coefficient, C = Ci/Cd where Ci is the interfacial drag
% tan_theta :   (scalar) tan of bed slope
% lambda    :   (scalar) dimensionless ice velocity lambda = (u_ice * Cd *L) / (u_infty * St * tau * c). = 1/M from the paper
% u_tidal   :   (scalar) dimensionless tidal contribution to the boundary layer velocity
% timescale :   (scalar) timescale of the problem (used to calculate the current tidal contribution), in days
% N 	    :   (integer) number of grid points inside intrusion length approx
% M 	    :   (integer) number of approx intrusion lengths
% dt        :   (scalar)  timestep
% t_end     :   (scalar)  end time
% t_store   :   (scalar)  time step at which solution is store#
% t_plot    :   (scalar)  time per plot

%% Outputs
% x 	    :    (1 x n) scalar: spatial grid, n: number of grid points
% t_out     :    (m x 1) scalar: time output points, m: number of time output points
% h_out     :    (m x n) scalar: output channel thickness at time output point
% h1_out    :    (m x n) scalar: output thickness of fresh layer at time output points
% u1_out    :    (m x n) scalar: output velocity of fresh layer at time output points
% Fr_out    :    (m x n) scalar: output Froude number at time output points
% melt_out  :    (m x n) scalar: output melt rate at time output points
% dhdx_out  :    (m x n) scalar: output channel slope
% dhdt_out  :    (m x n) scalar: output channel time derivative
% idxmax    :    integer       : highest index which sol touches

%%
%% Setup
%%

%comput quantities
n        = N*M;         %total number of grid points
lp       = -(1 - 1/ 4 / Fr_infty^2 - 3/4 * Fr_infty^(2/3)); %approximation to intrusion length for no interfacial drag, flat bed and flat channel
dx       = lp/N;
x        = -(dx:dx:n*dx); %grid points
n_timesteps     = floor(t_end/dt);   %number of timesteps
n_steps_plot    = floor(t_plot/dt);  %number of timesteps per plot
%n_steps_store   = floor(t_store/dt); %number of timesteps per storage


% Initial conditions
t_now = 0;
h_now = ones(1,length(x));

%storage
%sz       = 1 + floor(n_timesteps/n_steps_store); %number of storage points (also store first timestep)
sz       = length(idx_out);
h_out    = zeros(sz, length(x)); %channel thickness
h1_out   = zeros(size(h_out));   %layer 1 thickness
u1_out   = zeros(size(h_out));   %layer 1 velocity
Fr_out   = zeros(size(h_out));   %Froude number
melt_out = zeros(size(h_out));   %melt rate
dhdx_out = zeros(size(h_out));   %solution slope
dhdt_out = zeros(size(h_out));
t_out    = zeros(sz,1);

%%
%% Timestepping
%%
count = 1;
idxmax = 1; %maximum index of x which is included in the solution (used for selecting appropriate axes when plotting)
for i = 1:n_timesteps+1
        %% compute the slope of the channel
        h_now = (h_now); 
        dhdx_now = compute_slope(h_now, dx);

        %% create splines of the channel width and slope (ODE solver needs slope and width at finer resolution)
        h_spline = spline(x,h_now);
        dhdx_spline = spline(x,dhdx_now);

        %% solve the two-layer equation for h1
        odefun  = @(x,h1) two_layer_RHS(x,h1,C,S, Fr_infty, h_spline, dhdx_spline);
        options = odeset('RelTol', 1e-3,...
                         'AbsTol', 1e-4,...
                        'Events', @(x,h1) evolve_shape_event_function(x,h1, h_spline));
        BC     = (Fr_infty*(1 + 1e-4))^(2/3) ;            %boundary condition: h1 at entry corresponds to Fr = 1
        sol    = ode15s(odefun, [0,min(x)], BC, options); %solve equations
       
        %% put the solution onto the grid
        idx = find(x > min(sol.x), 1,'last'); %x(1:idx) ARE in the solution grid (those points inside the wedge)
        idxmax = max([idxmax, idx]);
        x_solgrid = x(1:idx);        %x values which are in the solution grid
        h1_now = h_now;              %layer one equal to channel width...
        h1_now(1:idx) = deval(sol, x_solgrid); %...except in salt wedge
       
        %% compute associated quantities
        Fr_now    = Fr_infty .*  (h1_now).^(-3/2);
        u1_now    = 1./h1_now;
        
        %compute the tidal component
        t_now_dimensional = t_now * timescale;
        tidal_timescale_dimensional = 0.5; %tidal timescale in days
        tidal_timescale_dimensional_solar = 28; %tidal timescale in days

        u1_now    = u1_now + u_tidal*sin(2*pi*t_now_dimensional/tidal_timescale_dimensional)*sin(2*pi*t_now_dimensional/tidal_timescale_dimensional_solar);
        u1_now    = abs(u1_now); % bl velocity is positive

        %compute the melt
        phi_now   = h1_now./h_now;
        melt_now  = u1_now .* (1 - phi_now); %melt rate
        dhdt_now  = melt_now - lambda * (dhdx_now); %kinematic condition, no upwinding
        %dhdt_now  = melt_now - lambda * ((2 - 918/1028)*tan_theta / Cd + (918/1028) *dhdx_now);

     %% store the solution
        if any(i == idx_out)%mod(i-1, n_steps_store) == 0
                h_out(count,:)    = h_now;
                h1_out(count,:)   = h1_now;
                u1_out(count,:)   = u1_now;
                Fr_out(count,:)   = Fr_now;
                melt_out(count,:) = melt_now;
                t_out(count,:)    = t_now;
                dhdt_out(count,:) = dhdt_now;
                dhdx_out(count,:) = dhdx_now;
                count = count + 1;
                fprintf('solution stored at t = %.3f \n', t_now)
         end

%         %% plot the solution
%         if (mod(i-1,n_steps_plot)  == 0) && (t_plot <= t_end)
%                 figure(1); clf;
%                 plot_config(x, h_now, h1_now,u1_now,Fr_now, melt_now,t_now,dhdt_now,dhdx_now,...
%                            C, S,lambda)
%         %       pause
%                 drawnow
%         end
%  
        %% timestep the channel width
        h_prev = h_now;
        h_now1 = h_prev + dt * dhdt_now;
        h_now(1:end-1) = h_prev(1:end-1) - dt * lambda * (h_prev(1:end-1) - h_prev(2:end))/dx +  dt*melt_now(1:end-1); %upwinding scheme
        h_now(end)     = h_prev(end) - dt * lambda * (-1/2 *h_prev(end) + 2*h_prev(end-1) - 3/2 *h_prev(end-2))/dx +  dt*melt_now(end); %one sided fd for the final point
       
        t_now = t_now + dt;

        %% break if the intrusion is too long
        if idx > 0.9*length(x)
            break
        end
end




end





