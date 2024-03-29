% Generate the data used to make supplementary figure 10c of the manuscript.
%
% Produce a map of steady intrusion length as a function of F and dT for a
% given tidal velocity 

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
u_tidal= [1, 5, 10,50, 100]; %dimensionless tidal velocity 
xbig   = 1e5;                    
xeps   = 1e-4;                   

% 
% setup
%
intrusion_length = nan(length(dT), length(F), length(u_tidal));
nruns            = numel(intrusion_length);
count = 1;

% 
% do the solves
%
for iu = 1:length(u_tidal)
for idT = 1:length(dT)

    
    for iF = 1:length(F)
        if mod(count, 1000) == 0
            disp(['completed ', num2str(count), ' of ', num2str(nruns), ' solves'])
        end

        % get the steady solution
        [x2,y] = get_steady_problem_solution_consttidal(dT(idT),F(iF), C, S,u_tidal(iu),  xeps, xbig);

        if x2(end)==xbig
            L = nan;
            x = x2;
        else
            L = x2(end);
            x = x2-L;
        end 
        
        intrusion_length(idT, iF,iu) = L;
        count = count + 1;
    end
end
end


%
%% save the data
%
save('data-for-figures/supfigure10c_data.mat', 'dT', 'F', 'C', 'S','intrusion_length', 'u_tidal');