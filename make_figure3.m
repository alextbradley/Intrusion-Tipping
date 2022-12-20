% Make figure 3 in the ms, showing the flat steady intrusion length as a
% function of F and dT, with tipping point contours for different values of
% C. 
%
% 22/11/22 ATB, MIT license.
clear
addpath('functions')
%
% load the data
% 
f = load("data-for-figures/figure3_data.mat");
figure(1); clf; hold on; box on; grid on

%
% fill in the background
%
yl = [0.1, 10];
xl = [1e-1, 0.901]; %x and y limits of final plot
bgc = [67,85,90]/100; %background color
fill([xl, flip(xl)], [yl(1), yl(1), yl(end), yl(end)], bgc , 'FaceAlpha', 0.2)

%
% add map where finite intrusion exists
%
contourf(f.F, f.dT, log10(f.intrusion_length), 30, 'linestyle', 'none');
cc = colorbar; 
clabels = 10.^(-3:1:1);  
clabstring  =  ["10^{-3}", "10^{-2}", "10^{-1}", "10^{0}", "10^{1}"];
caxis(log10(clabels([1 end]))); 
set(cc,'Ticks',log10(clabels),'TickLabels', clabstring); 
cmap = cmocean('tempo');
colormap(cmap(20:end,:)); %remove final white shades
cc.Label.String = 'dimensionless intrusion length';
cc.FontSize = 15;


%
% add the tipping point curves for different C
%
colmap = cmocean('matter',6);
colmap = colmap((length(colmap)-3):end, :);
colmap(1,:) = ([0.9490    0.5855    0.4038] + [ 0.8874    0.3829    0.3249])/2;
for i = 1:3
    plot(f.FF, f.dTc(i,:),  'color', colmap(i,:), 'linewidth', 2.2);
end

%
% add text labels
%
txt1 = text(0.11,8.5, 'unbounded intrusion', 'FontSize', 16, 'FontName', 'Helvetica', 'Interpreter', 'none');
txt2 = text(0.63,0.12, 'bounded intrusion', 'FontSize', 16, 'FontName', 'Helvetica', 'Interpreter', 'none');
txt3 = text(0.565, 5.1, '$C = 0.01$', 'Interpreter', 'Latex', 'Rotation', 45, 'Color', colmap(3,:), 'FontSize', 16);
txt4 = text(0.7, 5.3, '$C = 0.1$', 'Interpreter', 'Latex', 'Rotation', 48, 'Color', colmap(2,:), 'FontSize', 16);
txt5 = text(0.84, 5.9, '$C = 1$', 'Interpreter', 'Latex', 'Rotation', 55, 'Color', colmap(1,:), 'FontSize', 16);




%
% tidy things
%
set(gca, 'YScale', 'log')
ylim(yl)
xlim(xl)
xlabel('$F$', 'FontSize', 16, 'Interpreter','latex')
ylabel('$M$', 'FontSize', 16, 'Interpreter','latex')
ax = gca; ax.XTick = 0.1:0.2:0.9;
ax.FontSize = 16;
fig = gcf;
fig.Position(3:4) = [560,420];
