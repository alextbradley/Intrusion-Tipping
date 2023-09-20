% Make supplementary figure 6, showing histograms of grounding line velocs 
% for the selected ice shelves

%% Preliminaries
%clear
addpath('functions');

%% Loop over shelves
figure(1); clf;
fnames = ["Amery", "Filchner", "Getz", "Larsen", "PIG", "PopeSmithKohler", "Ross", "Ronne", "Thwaites"];
titles =  ["Amery", "Filchner", "Getz", "Larsen", "Pine Island", "Pope, Smith, Kohler", "Ross", "Ronne", "Thwaites"];
for iF = 1:length(fnames)
    fshelf = load(strcat('./data-for-figures/shelves/', fnames(iF), '.mat'));

    subplot(3,3,iF);
    idx = ~isnan(fshelf.velocs);
    histogram(log10(fshelf.velocs(idx)),'Normalization',  'probability', 'FaceColor', [0, 33, 204]/256)
    hold on
    %mean(fshelf.velocs(idx))
    plot( log10(mean(fshelf.velocs(idx)))*[1,1],[0,1], 'k--', 'linewidth', 1.25)
    title(titles(iF), 'FontName', 'Arial');
    xlim(log10([10,5000]))
    ylim([0,0.5])
    %set(gca, 'XScale', 'log')
    xlabel('velocity (m/yr)')
    ylabel('density')
    ax = gca;
    ax.FontSize = 12;
    ax.XTick = [1,2,3];
    ax.XTickLabels = {'10^1', '10^2', '10^3'};

end %end loop over fnames

fig = gcf;
fig.Position(3:4) = [800, 800];
