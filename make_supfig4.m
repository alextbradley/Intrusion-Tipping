% Make supplementary figure 4, showing ice velocities and grounding
% line/ice front points for each of the ice shelves considered.

%% Preliminaries
%clear
addpath('functions');

%% load the big data file
f = load('data-for-figures/all-Ant-data.mat');

%% Loop over files
figure(1); clf;
fnames = ["Amery", "Filchner", "Getz", "Larsen", "PIG", "PopeSmithKohler", "Ross", "Ronne", "Thwaites"];

for iF = 1:length(fnames)
    fshelf = load(strcat('./data-for-figures/shelves/', fnames(iF), '.mat'));

    in = fshelf.in;
    [rId, cId] = find( in ) ; %rows and columns of non-zero entries
    xminidx = min(rId);
    xmaxidx = max(rId);
    yminidx = min(cId);
    ymaxidx = max(cId); %indices of min and max in rows and columns

    %extend these to make them square 
    xl = xmaxidx - xminidx;
    yl = ymaxidx - yminidx;
    if yl >= xl
        xminidx = xminidx - floor((yl - xl)/2);
        xmaxidx = xminidx + yl;
    else
        yminidx = yminidx - floor((xl - yl)/2);
        ymaxidx = yminidx + xl;

    end

    xs = f.x(xminidx:xmaxidx);
    ys = f.y(yminidx:ymaxidx);
    vs = f.vx(xminidx:xmaxidx, yminidx:ymaxidx);
    us = f.vy(xminidx:xmaxidx, yminidx:ymaxidx);
    tfmax = f.tf_max(xminidx:xmaxidx, yminidx:ymaxidx);
    speeds   = sqrt(us.^2 + vs.^2);                            %ice velocities

    subplot(3,3,iF); hold on
    contourf(xs, ys, speeds', 10, 'linestyle', 'none');
    clim([0, 4000])
    axis equal
    %xticks([]);
    %yticks([]);

    % add the gl and front points
    in = fshelf.in;
    [rId, cId] = find( in ) ; %rows and columns of non-zero entries
    xminidx = min(rId);
    xmaxidx = max(rId);
    yminidx = min(cId);
    ymaxidx = max(cId); %indices of min and max in rows and columns (back to what we actually use in the calculation)
    isgl = f.isgl(xminidx:xmaxidx,yminidx:ymaxidx);
    isfront = f.isfront(xminidx:xmaxidx,yminidx:ymaxidx);
    xs = f.x(xminidx:xmaxidx);
    ys = f.y(yminidx:ymaxidx);
    [xxs,yys] = meshgrid(xs,ys);
    xxs = xxs';
    yys = yys';
    idxgl = (isgl == 1);
    idxfr = (isfront == 1);

    scatter(xxs(idxgl),yys(idxgl), ones(length(xxs(idxgl))), 'k')
    scatter(xxs(idxfr),yys(idxfr), ones(length(xxs(idxfr))), 'r')
    drawnow ; pause
end %end loop over fnames

fig = gcf;
fig.Position(3:4) = [800, 800];
