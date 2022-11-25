function [x,y] = get_steady_problem_solution(dT, F, C, S, xeps, xbig)
% Return the solution to the steady state problem on the domain [xeps,
% xbig]

yeps = [1+2/3.*dT*(2*F^2*C/(1-F^2))^(1/2)*xeps^(3/2);
    1-(2*F^2*C/(1-F^2))^(1/2)*xeps^(1/2)]; %initial condition at x = xeps

[x,y] = ode45(@(x,y) ode_fun_steady(x,y,dT, F, C, S),[xeps xbig],yeps,...
    odeset('Events',@(x,y) events_fun_steady(x,y,F),'abstol',1e-7,'reltol',1e-7));
