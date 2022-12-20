% Generate the data used to make figure 3 of the manuscript.
%
% Produce a map of steady intrusion length as a function of F and dT.

%
% Preliminaries
%
clear
addpath('functions/');

%
% Parameters
%
N      = 300;  %number of points in both directions
lambda = 10.^linspace(-1,1,N); %inverse dT values
dT     = 1./lambda;
F      = linspace(0.05,0.9,N);     % F values
C      = 0.1;                    % drag cofficient
S      = 0;                      % dimensionless slope
xbig   = 1e5;                    
xeps   = 1e-4;                   

% 
% setup
%
intrusion_length = nan(length(dT), length(F));
nruns            = numel(intrusion_length);
count = 1;

% 
% do the solves
%
for idT = 1:length(dT)
    for iF = 1:length(F)
        if mod(count, 1000) == 0
            disp(['completed ', num2str(count), ' of ', num2str(nruns), ' solves'])
        end

        % get the steady solution
        [x2,y] = get_steady_problem_solution(dT(idT),F(iF), C, S, xeps, xbig);

        if x2(end)==xbig
            L = nan;
            x = x2;
        else
            L = x2(end);
            x = x2-L;
        end 
        
        intrusion_length(idT, iF) = L;
        count = count + 1;
    end
end

%
%% get the critical dT as a function of L for different C
%
FF  = linspace(0.05, 0.9, N);
CC  = [0.01, 0.1, 1];
dTc = nan(length(CC), length(FF));
count = 1;
for iF = 1:length(FF)
    for iC = 1:length(CC)
        dTc(iC,iF) = get_critical_M(S,CC(iC),FF(iF));
        if mod(count, 100) == 0
            disp(['completed ', num2str(count), ' of ', num2str(numel(dTc)), ' solves'])
        end
        count = count + 1;
    end
end

%
%% save the data
%
save('data-for-figures/figure3_data.mat', 'dT', 'F', 'C', 'S','intrusion_length', 'FF', 'CC', 'dTc');