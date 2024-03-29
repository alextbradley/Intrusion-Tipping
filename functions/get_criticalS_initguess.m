function Sc = get_criticalS_initguess(dT,C,F,guess)
% Determine the critical S for continuous intrusion as a function of the
% dimensionless temperature difference dT, far field Froude number F, and
% dimensionless drag coefficient C, with initial guess provided. 

%
% parameters
%
nitermax = 20; %maximum number of iterations
tol = 1e-3;    %solution tolerance
xeps = 1e-3; %initial condition location
xbig = 1e3;  %domain length

%
% initial guesses
%
lb = guess - .1;
ub = guess + .1;

%
% check that lb and ub are actually upper bounds
%
is_ub = 0;
while ~(is_ub)
    %get the solution to steady problem
    [x,~] = get_steady_problem_solution(dT, F, C, ub, xeps, xbig);

    if  x(end)  == xbig %we do have cts intrusion, so this is a legit upper bound
        is_ub = 1;
    else
        % disp('Initial upper bound guess was not an upper bound, extending search');
        ub = ub + 1;
    end
    ub;
end

is_lb = 0;
while ~(is_lb)
    %get the solution to steady problem
    [x,~] = get_steady_problem_solution(dT, F, C, lb, xeps, xbig);

    if  x(end)  < xbig %we don't have cts intrusion, so this is a legit lower bound
        is_lb = 1;
    else
        %disp('Initial lower bound guess was not an lower bound, extending search');
        lb = lb - 1;
    end
    lb;
end

%
% perform the bisection
%
niter = 1;   %seed 
while (abs(ub - lb) > tol && (niter < nitermax))
    S_guess = (ub + lb)/2;

    %solve equations with this value of dT
    [x,~] = get_steady_problem_solution(dT, F, C, S_guess, xeps, xbig);

    %update guess
    if  x(end)  == xbig %continuous intrusion, so we have a new upper bound
        ub = S_guess;
    else %finite intrusion, so new lower bound
        lb = S_guess;
    end
    niter = niter + 1;

end
Sc = S_guess;
end

