% Make supplementary figure 10 of the manuscript, showing the intrusion
% length as a function of time for different values of the tidal velocity.
%
% 06/11/23, ATB (aleey@bas.ac.uk). MIT licence
%
%% Preliminaries
addpath('functions');

fig = figure(1); clf; hold on; box on;
fig.Position(3:4) = [1200, 460];
fs = 14; %fontsize

positions = [0.08, 0.75, 0.4, 0.22;
            0.08, 0.1,   0.4, 0.55;
            0.55, 0.1,   0.4, 0.87];

for i = 1:3
    ax(i) = subplot('Position', positions(i,:));
    hold(ax(i), 'on');
    box(ax(i), 'on');
end

%% Panel a tidal forcing

t = linspace(0,28,1e3); %days

u_tidal = 0.02;
    tidal =  sin(2*pi*t / 0.5).*sin(2*pi*t/28);

plot(ax(1), t,tidal, 'linewidth', 1);
ax(1).XLabel.String = 'time (days)';
ax(1).YLabel.String = 'dimensionless tidal velocity';
ax(1).YLim = [-1.2,1.2];
ax(1).XLim = [0,28];

%% Panel b: intrusion distance as a function of time for different tidal velocities
% figure setup
basecolor = [0, 0,1; %cold
    1,0, 0]; %warm

colmap = nan(2,5,3);

% load the data
load('data-for-figures/supfigure10b_data.mat', 'intrusion_lengths', 'times', 'intrusion_length_consttidal', 't_consttidal');
sz = size(intrusion_lengths);
u_inf = 0.01;  %m/s, far field velocity
u_tidal_dimensional = [0.001,0.0025, 0.005, 0.01,0.025, 0.05, 0.1]; %m/s 
u_tidal   = u_tidal_dimensional/u_inf; % dimensionless tidal velocity

cols = linspecer(sz(2));
cols = parula( sz(2)+3);
cols = cols(3:end-1,:);
%cols = lines(sz(2));
for i = 1:sz(1)
    for j = 1:sz(2)
        %colmap(i,j,:) = basecolor(i,:)*(1-j/10);
        colmap(i,j,:) = cols(j,:);

    end
end

%parameters
lengthscale = 25; %x lengthscale

%% Make the plot
for i = 1:sz(1)
    for j = 1:sz(2)

        t = cell2mat(times(i,j));
        intrusion_length = cell2mat(intrusion_lengths(i,j));

        %have to remove zeros padded at the end of time
        intrusion_length = intrusion_length(t>0);
        t = t(t>0);
        if i==2
            plot(ax(2), t, intrusion_length*lengthscale, 'linewidth', 1.75, 'Color',colmap(i,j,:));
            legendinfo{j} = ['$U_T/U_{\infty} = ' num2str(u_tidal(j)) '$'];
        else
            plot(ax(2), t, intrusion_length*lengthscale, '--', 'linewidth', 1.75, 'Color',colmap(i,j,:), 'HandleVisibility','off');

        end

%    drawnow
%         pause


    end
end

% add the figure 2 results
fig2data = load('data-for-figures/figure2_data.mat');
fig2sol = fig2data.sol;


for i = 1:2

        intrusion_length = nan(1,length(fig2sol(i).t));
        t = nan(1,length(fig2sol(i).t));


        for iT = 1:(length(fig2sol(i).t)-1)
            %if we're in second case, add first row in low alpha

            idx = find(fig2sol(i).h(iT,:)-fig2sol(i).h1(iT,:) == 0, 1, 'First');
            intrusion_length(iT) = -fig2sol(i).x(idx)*lengthscale;
            t(iT) = fig2sol(i).t(iT);
        end
        if i==2
            plot(ax(2), t(1:end-1), intrusion_length(1:end-1), 'linewidth', 2, 'Color','k');
        else
            plot(ax(2), t(1:end-1), intrusion_length(1:end-1), '--', 'linewidth', 2, 'Color','k');

        end
        %       


end

% % add the constant forcing results (sol_const)
intrusion_length_consttidal = intrusion_length_consttidal(t_consttidal>0);
t_consttidal = t_consttidal(t_consttidal>0);
plot(ax(2), t_consttidal, intrusion_length_consttidal*lengthscale, 'r', 'linewidth', 1.5)


ax(2).XLabel.String = 'time (days)';
ax(2).YLabel.String = 'intrusion distance (m)';
set(ax(2), 'YScale', 'log');
ax(2).XLim = [0,50]; 
legend(ax(2), legendinfo, 'Interpreter','latex', 'location', 'northeast');
ax(2).YLim = [0, 7500];

%% Panel c: final intrusion length regime diagram with contours for different tidal velocity values
f = load('data-for-figures/supfigure10c_data.mat', 'dT', 'F', 'C', 'S','intrusion_length', 'u_tidal');
idx = find(f.u_tidal == 10);
%
% fill in the background
%
yl = [0.1, 10];
xl = [1e-1, 0.901]; %x and y limits of final plot
bgc = [67,85,90]/100; %background color
fill(ax(3),[xl, flip(xl)], [yl(1), yl(1), yl(end), yl(end)], bgc , 'FaceAlpha', 0.2)

%
% add map where finite intrusion exists for u_tidal = 10;
%
intr = f.intrusion_length;
intr = squeeze(intr(:,:,idx));
%p = imagesc(f.F, f.dT, log10(intr));
%set(p, 'AlphaData',~isnan(intr));
contourf(f.F, f.dT, log10(intr),100, 'linestyle', 'none');
set(ax(3), 'YDir', 'normal');
cc = colorbar(ax(3)); 
clabels = 10.^(-3:1:1);  
clabstring  =  ["10^{-3}", "10^{-2}", "10^{-1}", "10^{0}", "10^{1}"];
caxis(log10(clabels([1 end]))); 
set(cc,'Ticks',log10(clabels),'TickLabels', clabstring); 
cmap = cmocean('tempo'); colormap(cmap(20:end,:)); %remove final white shades
%cmap = cmocean('ice'); colormap(cmap(50:end-20,:)); %remove final white shades
cc.Label.String = 'dimensionless intrusion length';
cc.FontSize = fs;


%
% add text labels
%
txt1 = text(ax(3),0.11,8.5, 'unbounded intrusion', 'FontSize', fs, 'FontName', 'Helvetica', 'Interpreter', 'none');
txt2 = text(ax(3), 0.59,0.12, 'bounded intrusion', 'FontSize', fs, 'FontName', 'Helvetica', 'Interpreter', 'none');

%
% add contours of critical intrusion
%
colmap= flipud(cmocean('amp', length(u_tidal)+1));
for i = 1:length(f.u_tidal)
    intr = f.intrusion_length;
    intr = squeeze(intr(:,:,i));
    intr(isnan(intr)) = -20;
    contour(f.F, f.dT, intr, [-10,-10], 'color', colmap(i,:), 'linewidth', 1.5)


end

%
% tidy things
%
set(ax(3), 'YScale', 'log')
ax(3).YLim = yl;
ax(3).XLim = xl;
xlabel('$F$', 'FontSize', fs, 'Interpreter','latex')
ylabel('$M$', 'FontSize', fs, 'Interpreter','latex')
ax(3).XTick = 0.1:0.2:0.9;

%% tidy
for i = 1:3
    box(ax(i), 'on');
    ax(i).FontSize = fs;
end

%% save
set(0, 'DefaultFigureRenderer', 'painters');
exportgraphics(gcf, 'supfig10.eps', 'ContentType', 'Vector')