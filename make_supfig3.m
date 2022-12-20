% Make supplementary figure 3, showing the bisection procedure for (a)
% critical M, (b) critical S, (c) critical F

% AB 21/11/22 MIT License
%% Preliminaries
clear
addpath('functions');

%% Get critical M data
S = 0;
C = 0.1;
F = 0.25; %parameters from figure 2 of main text
[Mc,ss] = get_criticalM(S,C,F);

%% plot
figure(1); clf;
subplot(1,3,1); hold on; box on
plot([0,12], F^(2/3)*[1,1], 'k--', 'linewidth', 1.5)
colmap = cmocean('balance',length(ss));
for i = 1:length(ss)

    p = plot(ss(i).x, ss(i).y(:,2), 'linewidth', 1.25);
    err = abs(ss(i).M  - Mc);
    errmax = max([abs(ss(1).M  - Mc), abs(ss(end).M  - Mc)]);
    alpha = 0.2 + (1-0.2)*(1-err/errmax);
    if ss(i).M < Mc
       
        p.Color = [0, 33/256, 203/256, alpha];
    else
       p.Color =[203/256, 0, 33/256, alpha];
    end

end
xlabel('$x-x_d$', 'Interpreter','latex')
ylabel('$h_1$', 'Interpreter','latex')
xlim([0,10])
ax(1)= gca; ax(1).FontSize = 14;
ax(1).YLim = [0.2,1.5];

%% Get critical S data
M = 0.3;
C = 0.1;
F = 0.25; 
[Sc,ss] = get_criticalS(M,C,F);

%% plot
subplot(1,3,2); hold on; box on
plot([0,20], F^(2/3)*[1,1], 'k--', 'linewidth', 1.5)
colmap = cmocean('balance',length(ss));
for i = 1:length(ss)

    p = plot(ss(i).x, ss(i).y(:,2), 'linewidth', 1.25);
    err = abs(ss(i).S  - Sc);
    alpha = 0.2 + (1-0.2)*(1-err/0.4);
    if ss(i).S < Sc
       
       p.Color = [0, 33/256, 203/256, alpha];
    else
       p.Color =[203/256, 0, 33/256, alpha];
    end

end
xlabel('$x-x_d$', 'Interpreter','latex')
ylabel('$h_1$', 'Interpreter','latex')
xlim([0,10])
ax(2)= gca; ax(2).FontSize = 14;
ax(2).YLim = [0.2,1.5];

%% Get critical F data
M = 0.3;
C = 0.1;
S = 0;  
[Fc,ss] = get_criticalF(M,C,S);

% plot
subplot(1,3,3); hold on; box on
plot([0,20], F^(2/3)*[1,1], 'k--', 'linewidth', 1.5)
colmap = cmocean('balance',length(ss));
for i = 1:length(ss)

    p = plot(ss(i).x, ss(i).y(:,2), 'linewidth', 1.25);
    err = abs(ss(i).F  - Fc);
   
    alpha = 0.2 + (1-0.2)*(1-err/0.8);
    if ss(i).F < Fc
       
       p.Color = [0, 33/256, 203/256, alpha];
    else
       p.Color =[203/256, 0, 33/256, alpha];
    end

end
xlabel('$x-x_d$', 'Interpreter','latex')
ylabel('$h_1$', 'Interpreter','latex')
xlim([0,10])
ax(3)= gca; ax(3).FontSize = 14;
ax(3).YLim = [0.2,1.5];

%% tidy
for i = 1:3
    ax(i).XLabel.FontSize = 16;
    ax(i).YLabel.FontSize = 16;
    ax(i).XLim = [0,12];
end