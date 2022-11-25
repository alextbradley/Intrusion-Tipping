function [pos, isterm, dir] = evolve_shape_event_function(x,h1, h_spline)
%Events function: stop the solver when h reaches 1
h = ppval(h_spline, x); %channel width here
pos = h1 - (h-1e-3); %flag for h1 --> h
isterm = 1; %set to 1 to terminate integration
dir = 0; %trigger if h increases through 0
end

