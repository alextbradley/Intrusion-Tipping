% make figure 4 of the manuscript showing
% (a) intrusion length as a function of slope S for F = 0.25, C = 0.1, and
% different values of dT
% (b) regime diagram of critical slope as a function of dT and F
% (c) regime diagram of critical Fr as a function of dT and S
%
%
% Preliminaries
% 
clear
addpath('functions/')

% 
% setup the figure
%
positions = [0.1, 0.78, 0.85, 0.2;
             0.1, 0.42, 0.85, 0.3;
             0.1, 0.05, 0.85, 0.3];

fig = figure(1); clf;
for i = 1:3
    ax(i) = subplot('Position', positions(i,:));
    box on ; hold on; 

end
fig.Units = 'pixels'; fig.Position(3:4) = [550, 750];


colmap_b = cmocean('balance');
colmap_b = colmap_b([1:128, 129:2:256],:); 

%%
%
% plot (a)
%
f = load("data-for-figures/figure4a_data.mat");


%
% fill in the background for pro/retro grade slopes and dividing line
%
yl = [1e-2, 10]; %plot y lims
xl = [min(f.S), 5];

fill(ax(1),[min(f.S), 0, 0, min(f.S)]*1.1, [4.2, 4.2, yl(2), yl(2)],colmap_b(40,:), 'LineStyle', 'none') 
fill(ax(1),[0,xl(2), xl(2),0]*1.1, [4.2, 4.2, yl(2), yl(2)],colmap_b(end-20,:), 'LineStyle', 'none') 
plot(ax(1),[0,0], yl, 'k--', 'linewidth', 1.5)

%
colmap_a = flip(cmocean('speed', length(f.dT)+2));
cmap = colmap_a;
for i = 1:(length(f.dT)-1)
    plot(ax(1), f.S, f.intrusion_length(i,:), 'Color',cmap(i+1,:), 'LineWidth', 2)
    idx = find(isnan(f.intrusion_length(i,:)), 1,'first') - 1; %find last non nan entry
    plot(ax(1), f.S(idx), f.intrusion_length(i,idx), 'ko', 'MarkerFaceColor', cmap(i+1,:)); %plot final point as circle - critical slope

    % plot dashed line down to 0
    plot(ax(1), [f.S(idx) , f.S(idx)], [1e-4, f.intrusion_length(i,idx)], '--','Color',cmap(i+1,:), 'LineWidth', 1.25 )
end


set(ax(1), 'YScale', 'log'); 
grid(ax(1), 'on')
ax(1).YLim = yl;
ax(1).XLim = xl;
ax(1).FontSize = 14; 
ax(1).XLabel.Interpreter = 'latex';
ax(1).XLabel.String = '$S$';
ax(1).XLabel.FontSize = 16;
ax(1).YLabel.Interpreter = 'latex';
ax(1).YLabel.String = '$L$';
ax(1).YLabel.FontSize = 16;

% sort out colorbar
c(1) = colorbar(ax(1));
colormap(c(1), cmocean('speed'));
c(1).Ticks = [0,0.5,1];
c(1).TickLabels = {"10^{-1}", "10^0", "10^1"};
c(1).Label.Interpreter = 'latex';
c(1).Label.String = '$\Delta T$';
c(1).Label.FontSize = 16;

% add the text of pro/retrograde
txtr = text(ax(1), -9.8, 7, 'prograde', 'interpreter', 'none', 'FontSize', 16, 'color', 'w');
txtb = text(ax(1), .5, 7, 'retrograde', 'interpreter', 'none', 'FontSize', 16,'color', 'w');

shg

%% 
%
% plot b
%
f = load("data-for-figures/figure4b_data.mat");
cc = f.critical_slope;
cc(cc<-10) = -10;
contourf(ax(2), f.F, f.dT, cc, 100, 'linestyle', 'none');
colormap(ax(2), colmap_b);
set(ax(2), 'YScale', 'log'); 
clim(ax(2),[-10,5])
c(2) = colorbar(ax(2));
c(2).Label.Interpreter = 'latex';
c(2).Label.String = '$S_c$';
c(2).Label.FontSize = 16;
ax(2).Position(3) = ax(1).Position(3);
c(2).Position(1) = c(1).Position(1);
c(2).Label.Position(1) = 2.5; 
ax(2).XLim = [min(f.F)-1e-3, max(f.F)+ 1e-4];
ax(2).YLim = [min(f.dT)-1e-4, max(f.dT) + 1e-1]; %slightly adjust x and y lims to make box visible
%c(2).Ticks = -10:5:10;
ax(2).FontSize = 14; 
ax(2).XLabel.Interpreter = 'latex';
ax(2).XLabel.String = '$F$';
ax(2).XLabel.FontSize = 16;
ax(2).YLabel.Interpreter = 'latex';
ax(2).YLabel.String = '$\Delta T$';
ax(2).YLabel.FontSize = 16;

% add the zero contour -- the line in figure 3
hold on
contour(ax(2), f.F, f.dT, f.critical_slope, [0,0], 'k', 'linewidth',1.25);

%%
% 
% plot c
f = load("data-for-figures/figure4c_data.mat");
contourf(ax(3), f.S, f.dT, f.critical_F', 100, 'linestyle', 'none');
c(3) = colorbar(ax(3));
c(3).Label.Interpreter = 'latex';
c(3).Label.String = '$F_c$';
c(3).Label.FontSize = 16;
ax(3).FontSize = 14; 
ax(3).XLabel.Interpreter = 'latex';
ax(3).XLabel.String = '$S$';
ax(3).XLabel.FontSize = 16;
ax(3).YLabel.Interpreter = 'latex';
ax(3).YLabel.String = '$\Delta T$';
ax(3).YLabel.FontSize = 16;

colmap_c = cmocean('matter');
colormap(ax(3), colmap_c(20:end,:));
clim([0, 0.9])
plot([0,0], [0.1, 10], 'w--', 'linewidth', 1.5)

%% add the ice shelf data
folder  = 'data-for-figures/shelves/';
shelves              = ["Amery", "Filchner", "Getz", "Larsen", "PIGfast", "PopeSmithKohler", "Ronne", "Ross", "Thwaites"];
names                = ["Amery", "Filchner", "Getz", "Larsen", "PIG", "PSK", "Ronne", "Ross", "Thwaites"];
horizontal_alignment = ["right", "left", "left","left","left","left","right","right","center"];
xshift               = [-0.25,-0.5, 0.25, 0.25,0.15, 0.25,-.25, 0.25,0];
yshift               = [0    ,0.20, 0   ,  0  , 0   ,  0,  0,     0.5,-0.1 ];

%setup storage
sz = size(shelves);
mean_tf_max = nan(sz);
mean_slope  = nan(sz);
mean_velocs = nan(sz);

%get mean data
for i = 1:length(shelves)
    fname = strcat(folder, shelves(i), '.mat');
    f_shelf = load(fname);

    %thermal forcing
    tf_max = f_shelf.tf_shelf;
    tf_max = tf_max(~isnan(tf_max));
    mean_tf_max(i) = mean(tf_max);

    %slope
    slope = f_shelf.slope;
    slope = slope(~isnan(slope));
    mean_slope(i)= mean(slope);

    %velocity
    velocs = f_shelf.velocs;
    velocs = velocs(~isnan(velocs));
    mean_velocs(i)= mean(velocs);

end

% compute dimensionless quantities
L = 335000; 
cc = 3974;
St = 5.9e-4;
Cd = 1e-2;
uinf = 0.01; % 1 cm/s
secs_per_year = 365*24*60^2; %ice velocities are in m/a
dT = mean_tf_max / (L/cc)  * St / Cd * uinf ./ mean_velocs * secs_per_year;
S = tan(mean_slope)/Cd;

%add points
for i = 1:max(sz)
    plot(ax(3), S(i), dT(i), 'ko', 'markerfacecolor', 'k',  'markeredgecolor', 'k', 'MarkerSize', 5)
    t(i) = text(S(i)+ xshift(i), dT(i)+yshift(i), names(i),'FontSize', 16, 'HorizontalAlignment',horizontal_alignment(i));%, 'VerticalAlignment', 'bottom');
end

% final tidying
set(ax(3), 'YScale', 'log')
ax(3).Position(3) = ax(1).Position(3);
c(3).Position(1) = c(1).Position(1);
c(3).Ticks = 0:0.2:0.8;
xlim([-10,5])
