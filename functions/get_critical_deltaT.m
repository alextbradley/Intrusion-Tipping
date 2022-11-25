function dTc = get_critical_deltaT(S,C,F)
% Determine the critical \Delta T for continuous intrusion as a function of
% the slope S, far field Froude number F, and dimensionless drag coefficient C

%
% parameters
%
nitermax = 20; %maximum number of iterations
tol = 1e-5;    %solution tolerance
niter = 1;
xeps = 1e-3; %initial condition location
xbig = 1e3;  %domain length

%
% initial guesses
%
lb = 1e-3;
ub = 10;

%
% check that lb and ub are actually upper bounds
%
is_ub = 0;
while ~(is_ub)
    %get the solution to steady problem
    [x,~] = get_steady_problem_solution(ub, F, C, S, xeps, xbig);

    if  x(end)  == xbig %we do have cts intrusion, so this is a legit upper bound
        is_ub = 1;
    else
        % disp('Initial upper bound guess was not an upper bound, extending search');
        ub = ub * 10;
    end
end

is_lb = 0;
while ~(is_lb)
    %get the solution to steady problem
    [x,~] = get_steady_problem_solution(lb, F, C, S, xeps, xbig);

    if  x(end)  < xbig %we don't have cts intrusion, so this is a legit lower bound
        is_lb = 1;
    else
        %disp('Initial lower bound guess was not an lower bound, extending search');
        lb = lb * 0.1;
    end
end

%
% perform the bisection
%
while (abs(ub - lb) > tol && (niter < nitermax))
    dT_guess = (ub + lb)/2;

    %solve equations with this value of dT
    [x,~] = get_steady_problem_solution(dT_guess, F, C, S, xeps, xbig);

    %update guess
    if  x(end)  == xbig %continuous intrusion, so we have a new upper bound
        ub = dT_guess;
    else %finite intrusion, so new lower bound
        lb = dT_guess;
    end
    niter = niter + 1;

end
dTc = dT_guess;
end

