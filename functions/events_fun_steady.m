function [value,isterminal,direction] = events_fun_steady(x,y,F)
% events function to detect when flow becomes critical for the steady
% solutions
    H = y(1);
    H1 = y(2);
    value = H1-F^(2/3);
    isterminal = 1;
    direction = 0;
end