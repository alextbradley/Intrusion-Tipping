function f= two_layer_RHS(x,h1,C,S, Fr_infty, h_spline, dhdx_spline)
% Return the right hand side of the two layer equation for a given x
h = ppval(h_spline, x);             %determine the channel width at this value of x
dhdx = ppval(dhdx_spline, x); %as above with channel slope
Fr = Fr_infty ./(h1^(3/2)); %Froude number
f = Fr^2 * (C * h/(h-h1) + 1) - (S + dhdx);
f = f/(Fr^2 - 1);
end



