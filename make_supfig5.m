% Make supplementary figure 5, showing histograms of the grounding line
% slope. 

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
    idx = ~isnan(fshelf.slope);
    histogram(fshelf.slope(idx),'Normalization',  'probability', 'FaceColor', [0, 33, 204]/256)
    hold on
    mean(fshelf.slope(idx))
    plot( mean(fshelf.slope(idx))*[1,1],[0,1], 'k--', 'linewidth', 1.25)
    title(titles(iF), 'FontName', 'Arial');
    xlim([-0.2, 0.2])
    ylim([0,0.4])
    xlabel('bed slope')
    ylabel('density')
    ax = gca;
    ax.FontSize = 12;

end %end loop over fnames

fig = gcf;
fig.Position(3:4) = [800, 800];
