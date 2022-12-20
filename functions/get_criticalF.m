function [Fc,ss] = get_criticalF(dT,C,S)
% Determine the critical F for continuous intrusion as a function of the
% dimensionless temperature difference dT, far field Froude number S, and
% dimensionless drag coefficient C

%
% parameters
%
Fc  = nan; %seed this
nitermax = 20; %maximum number of iterations
tol = 1e-4;    %solution tolerance
xeps = 1e-3; %initial condition location
xbig = 1e3;  %domain length

%
% initial guesses
%
lb = 0.1;
ub = 0.9; %F must be between these values

%
% if we have unbounded intrusion for F = ub, then return a nan
%

%get the solution to steady problem
[x,~] = get_steady_problem_solution(dT, ub, C, S, xeps, xbig);

if  x(end) == xbig %continuous intrusion
   Fc = ub;
   return %return the nan 
end


%
% if we have bounded intrusion for F = lb, then return a nan
%
% ub
% lb
% xeps
% xbig
% dT 
% C 
% S
[x,~] = get_steady_problem_solution(dT, lb, C, S, xeps, xbig);

if  x(end)  < xbig %we don't have cts intrusion, so this return
    Fc = lb; 
    return
end
lb;

%
% perform the bisection
%
ss = struct;
niter = 1;   %seed 
while (abs(ub - lb) > tol && (niter < nitermax))
    F_guess = (ub + lb)/2;

    %solve equations with this value of dT
    [x,y] = get_steady_problem_solution(dT, F_guess, C, S, xeps, xbig);
    ss(niter).x = x;
    ss(niter).y = y;
    ss(niter).F = F_guess;

    %update guess
    if  x(end)  == xbig %continuous intrusion, so we have a new lower bound
        lb = F_guess;
    else %finite intrusion, so new upper bound
        ub = F_guess;
    end
    niter = niter + 1;

end
Fc = F_guess;
end

