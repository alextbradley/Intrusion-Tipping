% Make supplementary figure 6, showing histograms of thermal forcing for
% the selected ice shelves

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
    idx = ~isnan(fshelf.tf_shelf);
    histogram(fshelf.tf_shelf(idx),'Normalization',  'probability', 'FaceColor', [0, 33, 204]/256)
    hold on
    mean(fshelf.tf_shelf(idx))
    plot( mean(fshelf.tf_shelf(idx))*[1,1],[0,1], 'k--', 'linewidth', 1.25)
    title(titles(iF), 'FontName', 'Arial');
    xlim([0,4])
    ylim([0,0.4])
    xlabel('thermal forcing (C)')
    ylabel('density')
    ax = gca;
    ax.FontSize = 12;

end %end loop over fnames

fig = gcf;
fig.Position(3:4) = [800, 800];
