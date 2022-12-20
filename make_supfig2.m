%Make supplementary figure 2 of the manuscript showing (a)--(b) the
%solution of the steady intrusion problem for (a0 M < M_c and (b) M > M_c,
%and (c)--(d) this solution in (h,h_1) space.
% AB 21/11/22 MIT License
%% Preliminaries
clear
addpath('functions');

%% Parameters
pp = struct;
pp(1).F    = 0.25;   %froude number on entry Fr_infty = Q_{infty} / sqrt(g' * H_{infty}^3)
pp(1).C    = 0.1;    %dimensionless interfacial drag
pp(1).S    = 0;	 	 %bed slope
pp(1).M    = 0.27;   %dimensionless ocean temperature (= M)
pp(1).xeps = 1e-3;
pp(1).xbig = 1e5;

%set two;
pp(2)      = pp(1);  %make a copy
pp(2).M    = 0.30; %adjust dT

%% Solve the problems
for i = 1:2
    [x,y] = solve_steady_system(pp(i));
    pp(i).x = x;
    pp(i).y = y;
end

%% Make plots
fig = figure(1); clf;
fig.Position(3:4) = [800, 600];
count = 1;
colmap = [0, 33, 203; 203, 0, 33]/256;
for i = 1:2
    ax(count) = subplot(2,2,i);
    hold on; box on
    plot(pp(i).x, pp(i).y(:,1), 'color', colmap(1,:), 'linewidth', 1.5)
    plot(pp(i).x, pp(i).y(:,2), 'color', colmap(2,:), 'linewidth', 1.5)
    plot(ax(count).XLim, pp(i).F^(2/3)*[1,1], 'k--', 'linewidth', 1.5)
    xlabel('$x-x_d$', 'Interpreter','latex')
    ylabel('$h, ~h_1$', 'Interpreter','latex')

    % add a legend in (a)
    if i ==1
        legend({'$h$', '$h_1$'}, 'interpreter', 'latex', 'location', 'northwest', 'FontSize', 14)
    else
        ax(count).YTick = 0:20:60;
    end
    count = count + 1;

    ax(count) =  subplot(2,2,i+2); hold on; box on
    plot(pp(i).y(:,1), pp(i).y(:,2), 'color', 'k', 'linewidth', 1.5)
    xl = ax(count).XLim;
    plot(xl, pp(i).F^(2/3)*[1,1], 'k--', 'linewidth', 1.5)
    ax(count).XLim = xl;
    xlabel('$h$', 'Interpreter','latex')
    ylabel('$h_1$', 'Interpreter','latex')
    count = count + 1;
end

% add zoom inset in (c),(d)
axnew = axes();
axnew.Position = [0.74, 0.64, 0.150, 0.14];
hold on; box on;
i =2;
plot(pp(i).x, pp(i).y(:,1), 'color', colmap(1,:), 'linewidth', 1.5)
plot(pp(i).x, pp(i).y(:,2), 'color', colmap(2,:), 'linewidth', 1.5)
plot(axnew.XLim, pp(i).F^(2/3)*[1,1], 'k--', 'linewidth', 1.5)
xlabel('$x-x_d$', 'Interpreter','latex')
ylabel('$h, ~h_1$', 'Interpreter','latex')
axnew.XLim = [0,10];
axnew.YLim = [0,2];
axnew.FontSize = 14;
axnew.XLabel.Position(2) = -0.1;
axnew.YTick = 0:2;
axnew.XTick = [0,10];

axnew2 = axes();
axnew2.Position = [0.74, 0.165, 0.150, 0.14];
hold on; box on;
xlabel('$h$', 'Interpreter','latex')
ylabel('$h_1$', 'Interpreter','latex')
plot(pp(i).y(:,1), pp(i).y(:,2), 'color', 'k', 'linewidth', 1.5)
axnew2.XLim = [1,5];
xl = axnew2.XLim;
plot(xl, pp(i).F^(2/3)*[1,1], 'k--', 'linewidth', 1.5)
axnew2.XLim = xl;
axnew2.FontSize = 14;
axnew2.XTick = [1,5];
axnew2.XLabel.Position(2) = -.2;

axnew2.YLim = [0,3];




for i = 1:4
    ax(i).FontSize = 14;
    ax(i).XLabel.FontSize = 16;
    ax(i).YLabel.FontSize = 16;

    if i < 3 
        ax(i).YLim = [0,1.5];
    end
    
end
ax(4).YTick = 0:20:60;
ax(4).YLim = [0,40];
ax(4).XLim = [0,60];
ax(4).XTick = 0:20:60;




function dydx = ode_fun(x,y,pp)
% steady-odes for channel thickness H and upper-layer thickness H_1
H = y(1);
H1 = y(2);
dHdx = pp.M*H1.^(-1).*(1-H1./H);
dH1dx = H1.^3./(pp.F^2-H1.^3).*( pp.F^2./H1.^3*(pp.C*H./(H-H1)+1)-pp.S-dHdx );
dydx = [dHdx; dH1dx];
end


function [value,isterminal,direction] = events_fun(x,y,pp)
% events function to detect when flow becomes critical
H = y(1);
H1 = y(2);
value = H1-pp.F^(2/3);
isterminal = 1;
direction = 0;
end

function [x,y] = solve_steady_system(pp)

yeps = [1+2/3./pp.M*(2*pp.F^2*pp.C/(1-pp.F^2))^(1/2)*pp.xeps^(3/2);
    1-(2*pp.F^2*pp.C/(1-pp.F^2))^(1/2)*pp.xeps^(1/2)]; % starting behaviour near singularity at H1 = H = 1
options = odeset('Events',@(x,y) events_fun(x,y,pp),'abstol',1e-8,'reltol',1e-8);
[x,y] = ode45(@(x,y) ode_fun(x,y,pp),[pp.xeps pp.xbig],yeps,options);

end