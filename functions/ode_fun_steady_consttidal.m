function dydx = ode_fun_steady_consttidal(x,y,dT, F, C, S,u_tidal)
% steady-odes for channel thickness H and upper-layer thickness H_1
    H = y(1);
    H1 = y(2);
    dHdx = dT*(0.4049*u_tidal + H1.^(-1)).*(1-H1./H);
    dH1dx = H1.^3./(F^2-H1.^3).*( F^2./H1.^3*(C*H./(H-H1)+1)-S-dHdx );
    dydx = [dHdx; dH1dx];
    %x
end
