% Make figure 2 in the manuscript, showing (a) intrusion length as a
% function of time, (b) snapshots of evolution of the channel and (c)
% thermal driving and melt rate

%
% Load the solution
%
%clear
load('data-for-figures/figure2_data.mat')

%% Preliminaries
addpath('functions');
close(figure(1))
set(0,'DefaultTextInterpreter','latex','DefaultAxesFontSize',14,'DefaultTextFontSize',14);

%% Plot prep
fig = figure(1); clf; fig.Position(3:4) = [1100,800];
fs = 16; %axis label fontsize

ncols = 2;    %number of columns
nsnapshots = 5; %number of snapshots of channel
nrows = 1 + nsnapshots + 2;    %number of rows

colgap = 0.07;%column gap
rowgap = 0.012;
width = 0.4;  %subplot width
startx = (1 - width*ncols - (ncols-1)*colgap)/2 ; %x start pt
starty = 1 - 0.2;%y start point (leave room for colorbar)
height = starty/(nrows+ 1);
%height = 0.135;
positions = zeros(4, nrows, ncols);
for p = 1:nrows
    for q = 1:ncols
        positions(:,p,q) = [startx + (q-1)*colgap + (q-1)*width, starty - (p-1)*height  - (p-1)*rowgap, width, height];

        %make the first row higher and bigger
        if p == 1
            positions(2,p,q) =  positions(2,p,q) +  0.05;
            positions(4,p,q) =  positions(4,p,q) +  0.03;
        end
        %make the final rows low and bigger(?)
        if p == nrows || p == (nrows - 1)
            positions(2,p,q) =  positions(2,p,q) -  0.025;
            %positions(4,p,q) =  positions(2,p,q) +  0.01;
        end

        if p == nrows
            positions(2,p,q) =  positions(2,p,q) -  0.01;
            %positions(4,p,q) =  positions(2,p,q) +  0.01;
        end



        %         % make the first row larger
        %         if p == 1
        %             positions(4,p,q) =  positions(4,p,q) + 0.15;
        %         end

        ax(p,q) = subplot('Position', squeeze(positions(:,p,q)));
        hold(ax(p,q), 'On'); box(ax(p,q), 'On');
    end
end
shg

%% Snapshots

% Find the indices of time points
tidx = zeros(1,nsnapshots);
iP = 1;
tt = sol(iP).t;
tout = [1,5,10,50,100];
for i = 1:length(tout)
    [~,idx] = min(abs(tt - tout(i)));
    tidx(i) = idx;
end

%
% Colormaps
%
ns = 50; %number shift in colormap
C = zeros(length(sol(iP).t)+ns, 3, 2);

%colormap for the finite intrusion
A = cmocean('ice', (length(sol(iP).t)+ns));
A = flipud(A);
C(:,:,1) = A;

%colormap for the finite intrusion
A = cmocean('solar', (length(sol(iP).t)+ns));
A = flipud(A);
C(:,:,2) = A;

alexred = [178,34,34]/255;
iceblue = [231, 255,254]/256;
basebrown = [209, 199,186]/256;

%
% channel shape evolutions
ylims = [20, 20;
    20, 30;
    20, 40;
    20, 80;
    20, 150];
xlims = [120,800];
xticks = [1,30,60,90, 120;
    1, 200,400, 600, 800];
ylimlo =-0.1*ylims;

for iP = 1:2
    lengthscale = 25; %dimensional horizontal legnthscale (m)
    zscale = 10;      %dimensional vertical legnthscale (cm)

    for iT = 1:length(tidx)
        i = tidx(iT);
        plot(ax(iT+1,iP), -sol(iP).x*lengthscale,zeros(size(sol(iP).x)), 'k', 'linewidth', 1.5) %base

        %fill the top in blue
        xx = -sol(iP).x*lengthscale;
        xxf = [xx,flip(xx)];
        yy = sol(iP).h(i, :)*zscale;
        yyf = [yy,ylims(iT,iP)*.999*ones(size(yy)) ];
        fill(ax(iT+1,iP),xxf,yyf, iceblue, 'LineStyle','none', 'FaceAlpha', 0.5);

        %fill the base in brown
        yyf = [ylimlo(iT,iP)*ones(size(yy)), zeros(size(yy))];
        fill(ax(iT+1,iP),xxf,yyf, basebrown, 'LineStyle','none', 'FaceAlpha', 0.5);


        %add top and interface
        plot(ax(iT+1,iP),-sol(iP).x*lengthscale,sol(iP).h(i, :)*zscale, 'linewidth', 1.75, 'color', C(i,:,iP))
        pl = plot(ax(iT+1,iP),-sol(iP).x*lengthscale,sol(iP).h(i,:)*zscale-sol(iP).h1(i,:)*zscale,'--', 'linewidth', 1.25, 'color', C(i,:,iP));
        pl.Color = [ C(i,:,iP), 0.25];


        ax(iT+1,iP).YLim = [ylimlo(iT,iP),ylims(iT,iP)];
        ax(iT+1,iP).XLim = [1,xlims(iP)];
        ax(iT+1,iP).XTick = xticks(iP,:);
        ax(iT+1,iP).XTickLabel = {};
        ax(iT+1,iP).YTick =[0,ylims(iT,iP)];
        ax(iT+1,iP).FontSize = 14;
        ax(iT+1,iP).YLabel.String = 'y (cm)';
        ax(iT+1,iP).YLabel.Interpreter = 'none';
        ax(iT+1,iP).YLabel.FontName = 'Helvetica';
    end
end

%
%% Intrusion distance evolution
%

% compute the intrusion distance
timescale = 1; %days
ylimst = xlims; %use from previous
for iP = 1:2
    for iT = 1:2:length(sol(iP).t)
        %if we're in second case, add first row in low alpha
        if iP == 2
            idx = find(sol(1).h(iT,:)-sol(1).h1(iT,:) == 0, 1, 'First');
            p = scatter(ax(1,iP),sol(1).t(iT)*timescale, -sol(1).x(idx)*lengthscale,  'MarkerFaceColor', C(iT,:,1), 'MarkerEdgeColor', 'k','LineWidth', 0.5);
            p.MarkerFaceAlpha = 0.3;
            p.MarkerEdgeAlpha = 0.3;
        end
    
        idx = find(sol(iP).h(iT,:)-sol(iP).h1(iT,:) == 0, 1, 'First');  
        scatter(ax(1,iP),sol(iP).t(iT)*timescale, -sol(iP).x(idx)*lengthscale,  'MarkerFaceColor', C(iT,:,iP), 'MarkerEdgeColor', 'k','LineWidth', 0.5);
    
    end

    % Add the steady state solution in first column
    if iP == 1
        L = -sol(iP).x(idx)*lengthscale;
        plot(ax(1,iP),sol(iP).t*timescale,  L*ones(size(sol(iP).t)), '--', 'Color', 'k', 'LineWidth',1.5);
    end

    ax(1,iP).XLim = [0, sol(iP).t_end*timescale];
    ax(1,iP).XLabel.String = 'time (days)';
    ax(1,iP).XLabel.Interpreter = 'none';
    ax(1,iP).XLabel.FontName = 'Helvetica';
    ax(1,iP).YLabel.String = {'intrusion','length (m)'};
    ax(1,iP).YLabel.Interpreter = 'none';
    ax(1,iP).YLabel.FontName = 'Helvetica';
    %ax(1,iP).YLabel.String = '$x_{\mathrm{int}}$ (m)';
    ax(1,iP).FontSize = 14;
    ax(1,iP).YLim = 1.1*ax(1,iP).YLim;
    ax(1,iP).XScale = 'log';
    ax(1,iP).YLim = [0, ylimst(iP)];
    ax(1,iP).XLim = [1e-1, 100];
    ax(1,iP).YTick = xticks(iP,1:2:end);
    ax(1,iP).YTickLabels{1} = '0';
end

% Add the steady state to the snapshots
for iP = 1:2
    for iT = 1:length(tidx)
        if iP == 1
            plot(ax(iT+1,1),[L,L], [0,100], 'k--', 'linewidth', 1.5)
        end

        %compute the intrusion distance (where h - h1 = 0 first)
        i = tidx(iT);
        idx = find(sol(iP).h(i,:)-sol(iP).h1(i,:) == 0, 1, 'First');
        plot(ax(iT+1,iP), -sol(iP).x(idx)*lengthscale,0, 'o', 'MarkerFaceColor', C(i,:,iP), 'MarkerEdgeColor', 'k','LineWidth', 0.5, 'MarkerSize',7)


    end
end

%% Evolution of velocity and melt
for iP = 1:2
    phi = sol(iP).h1./sol(iP).h;
    %come up with some values of temp
    T1 = -0.5; %local freezing temp
    if iP == 2 %warm case
        T2 = 2;
    else
        T2 = 2*sol(1).pp.dT / sol(2).pp.dT;
    end
    Tbar = (T1 * phi + T2 * (1- phi));
    T2 - T1
    for iT = 1:length(tidx)
        i = tidx(iT);
        plot(ax(nsnapshots+2,iP), -sol(iP).x*lengthscale,Tbar(i,:)-T1, 'linewidth', 1.5, 'color', C(i,:,iP))
        %plot(ax(nsnapshots+2,iP), -sol(iP).x*lengthscale,phi(i,:), 'linewidth', 1.5, 'color', C(i,:,iP))

        %add the intrusion as pt
        idx = find(sol(iP).h(i,:)-sol(iP).h1(i,:) == 0, 1, 'First');
        plot(ax(nsnapshots+2,iP), -sol(iP).x(idx)*lengthscale,Tbar(i,idx)-T1, 'o', 'MarkerFaceColor', C(i,:,iP), 'MarkerEdgeColor', 'k','LineWidth', 0.5, 'MarkerSize',7)

    end
    ax(nsnapshots+2,iP).XLim = [1,xlims(iP)];
    ax(nsnapshots+2,iP).XTick = xticks(iP,:);
    ax(nsnapshots+2,iP).YLabel.String = {'thermal';'driving (C)'};
    ax(nsnapshots+2,iP).YLabel.Interpreter = 'none';
    ax(nsnapshots+2,iP).YLabel.FontName = 'Helvetica';
    ax(nsnapshots+2,iP).XTickLabel = {};
    ax(nsnapshots+2,iP).FontSize = 14;
    ax(nsnapshots+2,iP).YLim = [0, 2.5];

end





% velocity
for iP = 1:2
    gprimed = 0.27;
    veloc_scale = sol(iP).pp.F*sqrt(gprimed*zscale);

    for iT = 1:length(tidx)
        i = tidx(iT);
        plot(ax(nsnapshots+3,iP), -sol(iP).x*lengthscale,sol(iP).u1(i, :)*veloc_scale, 'linewidth', 1.5, 'color', C(i,:,iP))
        %add the intrusion as pt
        idx = find(sol(iP).h(i,:)-sol(iP).h1(i,:) == 0, 1, 'First');
        plot(ax(nsnapshots+3,iP), -sol(iP).x(idx)*lengthscale,sol(iP).u1(i, idx)*veloc_scale, 'o', 'MarkerFaceColor', C(i,:,iP), 'MarkerEdgeColor', 'k','LineWidth', 0.5, 'MarkerSize',7)



    end



    % tidy
    ax(nsnapshots+3,iP).XLim = [1,xlims(iP)];
    ax(nsnapshots+3,iP).FontSize = 14;
    ax(nsnapshots+3,iP).YLim = [0.4,1];
    ax(nsnapshots+3,iP).XTick = xticks(iP,:);
    ax(nsnapshots+3,iP).XTickLabel{1} = '0';
    ax(nsnapshots+3,iP).XLabel.String = 'channel position (m)';
    ax(nsnapshots+3,iP).XLabel.Interpreter = 'none';
    ax(nsnapshots+3,iP).XLabel.FontName = 'Helvetica';
    ax(nsnapshots+3,iP).YLabel.String = {'velocity';'(m/s)'};
    ax(nsnapshots+3,iP).YLabel.Interpreter = 'none';
    ax(nsnapshots+3,iP).YLabel.FontName = 'Helvetica';

end


%% make colorbars
fig = figure(2);
clf;
c1 = colorbar;
colormap(c1, C(:,:,1));
c1.Label.String = 'time (days)';
c1.Label.Interpreter = 'none';
c1.Label.FontName = 'Helvetica';
c1.FontSize = 14;
c1.Label.FontSize = 16;
c1.Ticks = linspace(0, 1, 5);
c1.TickLabels = {'10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}'};
c1.Location = 'north';
axn = gca;
axn.Visible = 'off';
c1.Position(4) = 0.02;


axnew = axes;
anxew.Position = ax.Position;
c2 = colorbar;
colormap(c2, C(:,:,2));
c2.Ticks = [];
c2.Location = 'north';
axnew.Visible = 'off';
c2.Position(4) = 0.02;
c2.Position(2) = c1.Position(2)+ 0.02;
