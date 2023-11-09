% Make supplementary figure 8, showing the intrusion distance as function
% of time for different slopes.

load('data-for-figures/supfigure8_data_.mat');
sz = size(intrusion_lengths);
sz = sz(2);
lengthscale = 25; %x lengthscale
cd = 1e-2;
slopes_unscaled = [-0.05, -0.01,-0.005,-0.001, 0, 0.0008, 0.0016]; 
slopes_scaled = slopes_unscaled/cd;
colormap = cmocean('haline',length(slopes_scaled)+1);

figure(1);clf; hold on

for i = 1:sz
    intrusion_length = cell2mat(intrusion_lengths(i));
    t = cell2mat(times(i));
    plot(t, (intrusion_length), 'linewidth', 2, 'Color',colormap(i,:));

    legendinfo{i} = sprintf('slope = %.3f percent', slopes_unscaled(i)*100);
end

legend(legendinfo, 'interpreter', 'latex', 'location', 'northwest');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
box on
ax = gca;
ax.XLim = [0, 50];
ax.FontSize = 14;
ax.YLim = [5e-2, 1e2];
xlabel('time (days)', 'interpreter', 'none');
ylabel('intrusion length (m)', 'interpreter', 'none')
fig = gcf;
fig.Position(3:4) =[560,420];
