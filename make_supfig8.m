% Make supplementary figure 8, showing the critical slope as a function of
% F and dT for C = 0.1.
%
% 14/04/23, ATB (aleey@bas.ac.uk)
%
%
addpath('functions/')
figure(1); clf; hold on
data = load('data-for-figures/figure4b_data.mat');

p = imagesc(data.F,data.dT, data.critical_slope);

%contourf(data.F, data.dT, data.critical_slope,20 , 'linestyle', 'none');
ax = gca;
ax.CLim =[-5,5];
set(ax, 'YDir', 'normal');

c = colorbar;
c.Label.String = '$S_c$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 14;


xlim([0.1,0.9])
ylim([0.1, 10]);
set(gca, 'YScale', 'log')
box on
colormap(cmocean('delta'))
ax.FontSize = 12;
ax.XLabel.String = '$F$';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.String = '$M$';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.FontSize = 14;
ax.YLabel.FontSize = 14;
