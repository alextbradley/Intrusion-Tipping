function dydx = ode_fun_steady(x,y,dT, F, C, S)
% steady-odes for channel thickness H and upper-layer thickness H_1
    H = y(1);
    H1 = y(2);
    dHdx = dT*H1.^(-1).*(1-H1./H);
    dH1dx = H1.^3./(F^2-H1.^3).*( F^2./H1.^3*(C*H./(H-H1)+1)-S-dHdx );
    dydx = [dHdx; dH1dx];
    %x
end
