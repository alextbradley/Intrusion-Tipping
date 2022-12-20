% Make figure 2 in the manuscript. For dT = 0.29 and 0.3, show the (a)
% evolution of the channel width and interface, (b) boundary layer velocity
% (c) thermal driving as functions of space. Also plot in (d) the intrusion
% length as a function of time.

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
fig = figure(1); clf; fig.Position(3:4) = [1100,600];
fs = 16; %axis label fontsize

ncols = 2;    %number of columns
nrows = 4;    %number of rows

colgap = 0.07;%column gap
rowgap = 0.012;
width = 0.4;  %subplot width
startx = (1 - width*ncols - (ncols-1)*colgap)/2 ; %x start pt
starty = 1 - 0.4;%y start point (leave room for colorbar)
height = starty/(nrows);
height = 0.135;
positions = zeros(4, nrows, ncols);
for p = 1:nrows
    for q = 1:ncols
        positions(:,p,q) = [startx + (q-1)*colgap + (q-1)*width, starty - (p-1)*height  - (p-1)*rowgap, width, height];

        %make the final row lower
        if p == nrows
            positions(2,p,q) =  positions(2,p,q) - 0.08;
        end

        % make the first row larger
        if p == 1
            positions(4,p,q) =  positions(4,p,q) + 0.15;
        end

        ax(p,q) = subplot('Position', squeeze(positions(:,p,q)));
        hold(ax(p,q), 'On'); box(ax(p,q), 'On');
    end
end
shg

%% Make plot


ns = 30; %number shift in colormap
alexred = [178,34,34]/255;
for iP = 1:2

    tidx = 1:3:length(sol(iP).t);
    A = cmocean('ice', (length(sol(iP).t)+ns));
    A = A(1:length(sol(iP).t),:);
    %A   = fakeparula(length(sol(iP).t));
    %
    %first row: evolution of channel shape
    %
    plot(ax(1,iP), -sol(iP).x,zeros(size(sol(iP).x)), 'k', 'linewidth', 1.5) %base

    for i = tidx
        plot(ax(1,iP),-sol(iP).x,sol(iP).h(i, :), 'linewidth', 1, 'color', A(i,:))
        plot(ax(1,iP),-sol(iP).x,sol(iP).h(i,:)-sol(iP).h1(i,:),'--', 'linewidth', 1, 'color', A(i,:));

        %also work out the intrusion distance (where h - h1 = 0 first)
        idx = find(sol(iP).h(i,:)-sol(iP).h1(i,:) == 0, 1, 'First');
        sol(iP).xint(i) = sol(iP).x(idx);
        plot(ax(1,iP), -sol(iP).xint(i),0, 'o', 'MarkerFaceColor', A(i,:), 'MarkerEdgeColor', 'k','LineWidth', 0.5)

    end
    ax(1,iP).YLabel.String = '$z$';
    ax(1,iP).YLim = [-0.1,2];
    %ax(1,iP).XLim = [-10, 0];
    ax(1,iP).XTickLabel = {};
    ax(1,iP).FontSize = fs;


    %
    %colorbar
    %
    if iP == 1
        c(iP)= colorbar(ax(1,iP));
        c1 = c(iP);
        c1.Location = 'northoutside';
        c1.Colormap = A();
        c1.Position(2) = 0.91;
        c1.Position(1) = 0.5 - width/2;

        c1.Label.String = 'dimensionless time';
        c1.Label.Interpreter = 'none';
        c1.Label.FontName = 'Helvetica';
        c1.Label.Position(2) = 3.15;
        c1.FontSize = 14;
        c1.Label.FontSize = 16;
        c1.Ticks = linspace(0, 1, 5);
        c1.TickLabels = {'10^{-2}', '10^{-1}', '10^{0}', '10^{1}', '10^{2}'};
    end
    % c1.Ticks = [0, 1/2,1];
    % c1.TickLabels = num2str(sol(iP).t(end)/2), num2str(sol(iP).t(end))};


    %
    % Add the steady state solution in first column
    %
    if iP == 1
        [x2,y] = get_steady_problem_solution(sol(iP).pp.dT,sol(iP).pp.F, sol(iP).pp.C, sol(iP).pp.S, 1e-4, 1e4);
        H = y(:,1); %steady thickness
        H1 = y(:,2);
        L = x2(end); x = x2-L;
        %plot(ax(1,iP),x,H, 'linewidth', 1.5, 'Color',alexred)
        plot(ax(1,iP),-sol(iP).x,sol(iP).h(i, :), 'linewidth', 2, 'Color',alexred)

        plot(ax(1,iP),[-100,x(1)],[1,1], 'linewidth', 2,  'Color',alexred)
        plot(ax(1,iP),-sol(iP).x,sol(iP).h(i, :) - sol(iP).h1(i,:),'linewidth', 2,  'Color',alexred);
        %plot(ax(1,iP),x,H-H1,'linewidth', 1.5,  'Color',alexred);
        plot(ax(1,iP),[-x(1),10],[0,0], 'linewidth', 2 , 'Color',alexred)
    end
    %drawnow; pause

    %
    %row two: mean temperature in the channel
    %
    phi = sol(iP).h1./sol(iP).h;

    %come up with some values of temp
    T1 = -0.5; %local freezing temp
    if iP == 2 %warm case
        T2 = 2;
    else
        T2 = 2*sol(1).pp.dT / sol(2).pp.dT;
    end
    Tbar = (T1 * phi + T2 * (1- phi))/(T2 - T1);
    for i = tidx
        plot(ax(2,iP), -sol(iP).x,Tbar(i,:), 'linewidth', 1, 'color', A(i,:))
    end
    ax(2,iP).YLabel.String = '$T$';
    ax(2,iP).XTickLabel = {};
    ax(2,iP).YLim = [-0.5, 1];
    ax(2,iP).YTick = -0.5:0.5:1;
    txt{iP} = sprintf("ocean temp: %.1f C ", T2);


    %
    %row three: boundary layer velocity
    %

    for i = tidx
        plot(ax(3,iP), -sol(iP).x,sol(iP).u1(i, :), 'linewidth', 1, 'color', A(i,:))
    end
    ax(3,iP).YLabel.String = '$u^*$';
    ax(3,iP).XLabel.String = 'dimensionless channel position';
    ax(3,iP).XLabel.Interpreter = 'none';
    ax(3,iP).XLabel. FontName = 'Helvetica';
    ax(3,iP).YLim = [1,2.5];
    ax(3,iP).YTick = 1:2;
    %ax(3,iP).XTickLabel = {};


    %
    %row four: intrusion distance
    %
    %add the steady state distance
    if iP == 1
        pp = plot(ax(4,iP),sol(iP).t,  -sol(iP).xint(i)*ones(size(sol(iP).t)), '--', 'Color', alexred, 'LineWidth',1.4);
        %pp = plot(ax(4,iP),sol(iP).t, -L*ones(size(sol(iP).t)), '--', 'Color', alexred, 'LineWidth',.75);
    end
    for i = tidx
        plot(ax(4,iP),sol(iP).t(i), -sol(iP).xint(i),'o',  'MarkerFaceColor', A(i,:), 'MarkerEdgeColor', 'k','LineWidth', 0.5)

    end
    
    ax(4,iP).XLim = [0, sol(iP).t_end];
    ax(4,iP).XLabel.String = 'dimensionless time';
    ax(4,iP).XLabel.Interpreter = 'none';
    ax(4,iP).XLabel.FontName = 'Helvetica';
    ax(4,iP).YLabel.String = '$x_{\mathrm{int}}$';
    ax(4,iP).YLim = 1.1*ax(4,iP).YLim;
    ax(4,iP).XLabel.FontSize = fs;
    ax(4,iP).XScale = 'log';
    ax(4,iP).YLim(1) = 0;
    ax(4,iP).XLim = [1e-2, max(sol(iP).t)];

end


for p = 1:nrows
    for q = 1:ncols
        ax(p,q).FontSize = 12;
        ax(p,q).YLabel.FontSize = 14;
        ax(p,q).XLabel.FontSize = 14;
        if p < nrows
            ax(p,q).XLim = [0,10];
        end
        ax(p,q).FontSize = 14;
        ax(p,q).XLabel.FontSize = 16;
        ax(p,q).YLabel.FontSize = 16;
    end
end

%% add the text title
t1 = cell2mat(txt{1,1});
txt1 = text(ax(1,1), 0, 2.2, t1, 'FontSize', 18, 'Interpreter', 'none', 'FontName','Helvetica');
t2 = cell2mat(txt{1,2});
txt2 = text(ax(1,2), 10-3.4, 2.2, t2, 'FontSize', 18, 'Interpreter', 'none', 'FontName','Helvetica');

% %% make the correct size
% fig = gcf;
% fig.Units = 'centimeters';
% target_width = 18; %18 cm
% scale = fig.OuterPosition(3)/target_width;
% fig.Position(3:4) = fig.Position(3:4)*scale;


