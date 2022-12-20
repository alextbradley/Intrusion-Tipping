% Generate the data for figure 4 of the manuscript. This figure shows
% (a) intrusion length as a function of slope S for F = 0.25, C = 0.1, and
% different values of dT
% (b) regime diagram of critical slope as a function of dT and F
% (c) regime diagram of critical Fr as a function of dT and S

%
% Preliminaries
%
clear
addpath('functions')

%which parts to gendata for flags
gendata_partA = 0;
gendata_partB = 0;
gendata_partC = 1;


%% Figure part a
if gendata_partA
    %
    % parameters
    %
    N = 500;                    %number of slope points
    S = linspace(-10,1,N);   %slope values
    F = 0.5;                   %Froude number
    dT = [10, 5,2, 1, 0.5,0.2, 0.1, 0.001];  %dT values (consider final one to be the no melt case)
    C = 0.1;                    %drag coefficient

    xeps = 1e-3;
    xbig = 1e4; %steady problem solver limits

    %
    % get the critical slope values
    %
    critical_slope = nan(length(dT),1);
    for idT = 1:length(dT)
        critical_slope(idT) = get_criticalS(dT(idT),C,F);
    end

    %
    %generate the data
    %
    intrusion_length = nan(length(dT), N);
    for iS = 1:N
        for idT = 1:length(dT)
            % get the steady solution
            if S(iS) < critical_slope(idT)
                [x2,~] = get_steady_problem_solution(dT(idT),F, C, S(iS), xeps, xbig);

                if x2(end)==xbig
                    L = nan;
                else
                    L = x2(end);
                end
                intrusion_length(idT,iS) = L;
            end
        end
    end

    %
    % save the data
    %
    save('data-for-figures/figure4a_data.mat', 'dT', 'C', 'intrusion_length', 'S', 'F', 'critical_slope');
end
%% Part b: critical slope as a function of dT and F
if gendata_partB
    N      = 30;                  %number of output points
    C      = 0.1;                  %drag coefficient
    lambda = 10.^linspace(-1,1,N); %inverse dT values
    dT     = 1./lambda;
    F      = linspace(0.1,0.9,N); % F values

    count = 1;
    tic
    nruns = length(dT)*length(F);
    critical_slope = nan(length(F),length(dT));
    critical_slope_flat = nan(1,length(dT)*length(F));
    for idT = 1:length(dT)
        %idT
        for iF = 1:length(F)
            %iF

            if ((iF == 1) && (idT == 1)) %start bisection fresh
               critical_slope(iF,idT) = get_criticalS(dT(idT),C,F(iF)); 
            elseif iF == 1
                guess = critical_slope_flat(count - length(F));
                %pause
                critical_slope(iF,idT) = get_criticalS_initguess(dT(idT),C,F(iF),guess);
            else %use previous guess
                guess = critical_slope_flat(count - 1);
                %pause
                critical_slope(iF,idT) = get_criticalS_initguess(dT(idT),C,F(iF),guess);
            end

            %critical_slope(iF,idT) = get_criticalS(dT(idT),C,F(iF));
            critical_slope_flat(count) = critical_slope(iF, idT);

            %update progress
            if mod(count, 1) == 0
                disp(['completed ', num2str(count), ' of ', num2str(nruns), ' solves'])
            end
            count = count + 1;
        end
    end

    %
    % save the data
    %
    save('data-for-figures/figure4b_data.mat', 'dT', 'C', 'F', 'critical_slope');
toc
end


%% Part c: critical Fr as a function of dT and S
if gendata_partC
    N  = 50;     %number of output points
    S  = -10:0.25:5;
    dT = 10.^linspace(-1,1.4,N);
    C  = 0.1;

    nruns = length(dT)*length(S);
    critical_F = nan(length(S),length(dT));
    count = 1;
    for iS = 1:length(S)
        for idT = 1:length(dT)
        idT,iS
        
            critical_F(iS,idT) = get_criticalF(dT(idT),C,S(iS));

            %update progress
            if mod(count, 10) == 0
                disp(['completed ', num2str(count), ' of ', num2str(nruns), ' solves'])
            end
            count = count + 1;
        end
    end


    save('data-for-figures/figure4c_data.mat', 'dT', 'S', 'critical_F', 'C');


end 