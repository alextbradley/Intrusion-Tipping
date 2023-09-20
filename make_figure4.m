% make figure 4 of the manuscript showing
% (a) intrusion length as a function of slope S for F = 0.25, C = 0.1, and
%     different values of dT
% (b) regime diagram of critical slope as a function of dT and F
% (c) regime diagram of critical Fr as a function of dT and S
%
%
% Preliminaries
% 
clear
addpath('functions/')
set(0,'DefaultTextInterpreter','none')


% 
% setup the figure
%
positions = [0.12, 0.78, 0.73, 0.2;
             0.12, 0.51, 0.73, 0.2;
             0.12, 0.08, 0.73, 0.35];

fig = figure(1); clf;
for i = 1:3
    ax(i) = subplot('Position', positions(i,:));
    box on ; hold on; 

end
fig.Units = 'pixels'; fig.Position(3:4) = [550, 750];


colmap_b = cmocean('balance');
colmap_b = colmap_b([1:128, 129:2:256],:); 
ax(1).XLabel.Position(2) = 0.0035;
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
xl = [-7,3];

%
colmap_a = flip(cmocean('speed', length(f.dT)+2));
cmap = colmap_a;
for i = 1:(length(f.dT)-1)
    plot(ax(1), f.S, f.intrusion_length(i,:), 'Color',cmap(i+1,:), 'LineWidth', 2)
    idx = find(isnan(f.intrusion_length(i,:)), 1,'first') - 1; %find last non nan entry
    plot(ax(1), f.S(idx), f.intrusion_length(i,idx), 'ko', 'MarkerFaceColor', cmap(i+1,:)); %plot final point as circle - critical slope

    % plot dashed line upwards
    plot(ax(1), [f.S(idx) , f.S(idx)], [f.intrusion_length(i,idx), 10], '--','Color',cmap(i+1,:), 'LineWidth', 1.25 )
end


set(ax(1), 'YScale', 'log'); 
grid(ax(1), 'on')
ax(1).YLim = yl;
ax(1).XLim = xl;
ax(1).FontSize = 16; 
ax(1).XLabel.Interpreter = 'latex';
ax(1).XLabel.String = '$S$';
ax(1).XLabel.FontSize = 18;
ax(1).YLabel.Interpreter = 'latex';
ax(1).YLabel.String = '$L$';
ax(1).YLabel.FontSize = 18;
ax(1).YTick = [0.01, 0.1, 1,10];

% sort out colorbar
c(1) = colorbar(ax(1));
colormap(c(1), cmocean('speed'));
c(1).Ticks = [0,0.5,1];
c(1).TickLabels = {"10^{-1}", "10^0", "10^1"};
c(1).Label.Interpreter = 'latex';
c(1).Label.String = '$M$';
c(1).Label.FontSize = 18;

% add the text of pro/retrograde
fill(ax(1),[min(f.S), 0, 0, min(f.S)]*1.1, [4.2, 4.2, yl(2), yl(2)],colmap_b(40,:), 'LineStyle', 'none') 
fill(ax(1),[0,xl(2), xl(2),0]*1.1, [4.2, 4.2, yl(2), yl(2)],colmap_b(end-20,:), 'LineStyle', 'none') 
plot(ax(1),[0,0], yl, 'k--', 'linewidth', 1.5)
txtr = text(ax(1), -6.8, 7, 'prograde', 'interpreter', 'none', 'FontSize', 16, 'color', 'w');
txtb = text(ax(1), .5, 7, 'retrograde', 'interpreter', 'none', 'FontSize', 16,'color', 'w');

shg

%% 
%
% plot b
%
% % add the red and blue boxes
% fill(ax(2),[8.2, 8.2, yl(2), yl(2)],[-10, 0, 0, -10]*1.1,colmap_b(40,:), 'LineStyle', 'none') 
% fill(ax(2),[8.2, 8.2, yl(2), yl(2)],[0,3, 3,0]*1.1 , colmap_b(end-20,:), 'LineStyle', 'none') 
% 
% g = load("data-for-figures/figure4b_data.mat");
% 
% 
% idxF = round(linspace(1,30,5)); %indices corresponding to F = 0.1, 0.3, 0.5, 0.7, 0.9 (roughly)
% cmapbname = 'dense';
% colmap_b2 = cmocean('dense', 100);
% colmap_b2 = colmap_b2(20:end-20,:);
% colmapblines= colmap_b2(round(linspace(1,length(colmap_b2),5)),:);
% plot(ax(2), [0.01,10], [0,0],'k--','linewidth', 1.25)
% for iF = 1:5
%     ii = idxF(iF);
%     plot(ax(2), g.dT, g.critical_slope(ii,:), 'linewidth',1.5, 'color', colmapblines(iF,:));
% end
% 
% %add the data from previous figure
% idxF = 6;
% for i = 1:(length(f.dT)-1)
%     idx = find(isnan(f.intrusion_length(i,:)), 1,'first') - 1; %find last non nan entry
%     plot(ax(2), f.dT(i), f.S(idx), 'ko', 'MarkerFaceColor', cmap(i+1,:)); %plot final point as circle - critical slope
% 
% end
% 
% ax(2).XLabel.Interpreter = 'latex';
% ax(2).XLabel.String = '$M$';
% ax(2).XLabel.FontSize = 18;
% ax(2).YLabel.Interpreter = 'latex';
% ax(2).YLabel.String = '$S_c$';
% ax(2).YLabel.FontSize = 18;
% %ax(2).XLabel.Position(2) = 0.065;
% set(ax(2), 'XScale', 'log'); 
% ax(2).YLim = [-3,3];
% ax(2).XLim = [1e-1,10];
% c(2) = colorbar(ax(2));
% c(2).Label.Interpreter = 'latex';
% c(2).Label.String = '$F$';
% c(2).Label.FontSize = 18;
% c(2).Limits = [min(g.F), max(g.F)];
% c(2).Ticks = 0.1:0.2:0.9;
% ax(2).FontSize = 16; 
% colormap(c(2), colmap_b2);

f = load("data-for-figures/figure4c_data.mat");

dT_idx = round(linspace(1,length(f.dT)-9,7));
dTs = f.dT;
dTs = dTs(dT_idx); %dT values to plot
Ms = 1./dTs;
Fc = f.critical_F;
cmap = flip(colmap_a);
for iM = 1:length(Ms)
    plot(ax(2),f.S, Fc(:,dT_idx(iM)), 'color', cmap(iM+1,:), 'LineWidth', 1.25);
end
   
ax(2).XLabel.Interpreter = 'latex';
ax(2).XLabel.String = '$S$';
ax(2).XLabel.FontSize = 18;
ax(2).YLabel.Interpreter = 'latex';
ax(2).YLabel.String = '$F_c$';
ax(2).YLabel.FontSize = 18;
ax(2).FontSize = 16; 
ax(2).XLim = [xl(1)-0.02, xl(2)]; %make box visibile
grid(ax(2), 'on')
ax(2).YTick = 0:0.2:1;

%%
% 
% plot c
f = load("data-for-figures/figure4c_data.mat");
%contourf(ax(3), f.S, f.dT, f.critical_F', 100, 'linestyle', 'none');
imagesc(ax(3), f.S, f.dT, f.critical_F'); set(ax(3), 'YDir', 'Normal')
%pcolor(ax(3), f.S, f.dT, f.critical_F'); set(ax(3), 'YDir', 'Normal')
c(3) = colorbar(ax(3));
c(3).Label.Interpreter = 'latex';
c(3).Label.String = '$F_c$';
c(3).Label.FontSize = 18;
ax(3).FontSize = 16; 
ax(3).XLabel.Interpreter = 'latex';
ax(3).XLabel.String = '$S$';
ax(3).XLabel.FontSize = 18;
ax(3).YLabel.Interpreter = 'latex';
ax(3).YLabel.String = '$M$';
ax(3).YLabel.FontSize = 18;
%ax(3).XLabel.Position(2) = ax(2).XLabel.Position(2);
ax(3).XLim = [xl(1)-0.02, xl(2)]; %make box visibile
ax(3).YLim = [0.1, 20.1];

colmap_c = cmocean('matter');
colormap(ax(3), colmap_c(20:end,:));
caxis([0, 0.9])
plot([0,0], [0.1, 22], 'k--', 'linewidth', 1.5)

%% add a couple of contours corresponding to ice sheets
levs = [0.3,0.7];
for i = 1:2
    contour(ax(3), f.S, f.dT, f.critical_F', [levs(i), levs(i)] , 'linewidth', 1.5, 'linestyle', '--', 'linecolor', 0.75*[1,1,1]);
end
%% add the ice shelf data
folder  = 'data-for-figures/shelves/';
shelves              = ["Amery", "Filchner", "Getz", "Larsen", "PIGfast", "PopeSmithKohler", "Ronne", "Ross", "Thwaites"];
names                = ["Amery", "Filchner", "Getz", "Larsen", "PIG", "PSK", "Ronne", "Ross", "Thwaites"];
horizontal_alignment = ["right", "left", "left","left","left","left","right","right","center"];
xshift               = [-0.25,-0.5, 0.25, 0.25,0.15, 0.25,-.25, 0.25,0];
yshift               = [0    ,0.20, 0   ,  0  , 0   ,  0,  0,     0.5,-0.1 ];

% constants
L = 335000; 
cc = 3974;
St = 5.9e-4;
Cd = 1e-2; %1e-2
uinf = 0.01; % 1 cm/s 0.01
secs_per_year = 365*24*60^2; %ice velocities are in m/a

%setup storage
sz = size(shelves);
mean_tf_max = nan(sz);
mean_slope  = nan(sz);
mean_velocs = nan(sz);
std_tf_max = nan(sz);
std_slope  = nan(sz);
std_velocs = nan(sz);

%alternative method: work out S and M for all points, then take the mean
clear S_all M_all

%get mean data
for i = 1:length(shelves)
    fname = strcat(folder, shelves(i), '.mat');
    f_shelf = load(fname);

    %thermal forcing
    tf_max = f_shelf.tf_shelf;
    tf_max = tf_max(~isnan(tf_max));
    mean_tf_max(i) = median(tf_max);
    std_tf_max(i)  = std(tf_max);

    %slope
    slope = f_shelf.slope;
    slope = slope(~isnan(slope));
    mean_slope(i)= median(slope); 
    std_slope(i) = std(slope);
    uq_slope(i) = quantile(slope, [0.7]);
    lq_slope(i) = quantile(slope, [0.3]);
    

    %velocity
    velocs = f_shelf.velocs;
    velocs = velocs(~isnan(velocs));
    mean_velocs(i)= mean(velocs);
    mean_velocs(i)= median(velocs);
    invmean_velocs(i) = mean(1./velocs);
    std_velocs(i) = std(velocs);
    uq_velocs(i) = quantile(velocs, [0.7]);
    lq_velocs(i) = quantile(velocs, [0.3]);

    
    %compute all values of S and M
    S_all = tan(slope)/Cd;
    M_all = mean_tf_max(i) / (L/cc)  * St / Cd * uinf ./ velocs * secs_per_year;
    
%      dT(i) = mean(M_all);
%      S(i) = mean(S_all);
end

% compute dimensionless quantities
dT     = mean_tf_max / (L/cc)  * St / Cd * uinf ./ mean_velocs * secs_per_year;
dT_min = mean_tf_max / (L/cc)  * St / Cd * uinf ./ (uq_velocs) * secs_per_year;
dT_max = mean_tf_max / (L/cc)  * St / Cd * uinf ./ (lq_velocs) * secs_per_year;

S     = tan(mean_slope)/Cd;
S_min = tan(lq_slope)/Cd;
S_max = tan(uq_slope)/Cd;


%add points
for i = 1:max(sz)
    plot(ax(3), S(i), dT(i), 'ko', 'markerfacecolor', 'k',  'markeredgecolor', 'k', 'MarkerSize', 5)
    plot(ax(3), [S(i),S(i)], [dT_min(i), dT_max(i)], 'k', 'LineWidth', 0.5)
    
   plot(ax(3), [S_min(i),S_max(i)], [dT(i), dT(i)], 'k', 'LineWidth', 0.5)
    t(i) = text(S(i)+ xshift(i), dT(i)+yshift(i), names(i),'FontSize', 16, 'HorizontalAlignment',horizontal_alignment(i));%, 'VerticalAlignment', 'bottom');
   % drawnow
   % pause
end

% final tidying
set(ax(3), 'YScale', 'log')
ax(3).Position(3) = ax(1).Position(3);
c(3).Position(1) = c(1).Position(1);
c(3).Ticks = 0:0.2:0.8;
xlim(xl)
for i = [1,3]
    c(i).Position(1) = 0.86;
    ax(i).Position = positions(i,:); %make sure nothing has been upset
end
ax(1).XLabel.Position(1) = -2.5;
ax(3).XLabel.Position(1) = -2.5;
ax(2).XLabel.Position(1) = -2.5;

