%% Preliminaries
%data = load('data-for-figures/all-Ant-data.mat');
%speed = sqrt((data.vx).^2 + (data.vy).^2);
%u = data.vx;
%v = data.vy;
%set(0, 'DefaultFigureRenderer', 'painters');

front_clr = 0.5*[1,1,1];
gl_clr    = [0,0,0];           %colours for the front and grounding line

gl_clr = [230, 43, 43]/255;
front_clr = 0.1*[1,1,1];


shelves = ["Thwaites", "Ross", "Ronne", "PopeSmithKohler", "PIG", "Larsen", "Getz", "Filchner", "Amery"];
%figure(1); clf;
cmap = cmocean('tempo', 100);
cmap = cmap(1:end-20, :);


%%
for i = 9
    fig(i) = figure(i); clf;
    ax(i) = gca;
    fname = strcat("data-for-figures/shelves/", shelves(i) , ".mat");

    shelf = load(fname);

    [rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
    xminidx = min(rId); xmaxidx = max(rId);
    np = 40; % number of padding points
    xminidx = xminidx - np;
    xmaxidx = xmaxidx + np;
    yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
    yminidx = yminidx - np;
    ymaxidx = ymaxidx + np;


    dy = ymaxidx - yminidx;
    dx = xmaxidx - xminidx;
    dd = abs(dx - dy);
    if dy < dx %adjust so square
        yminidx = yminidx - floor(dd/2);
        ymaxidx = ymaxidx + floor(dd/2);
    else
        xminidx = xminidx - floor(dd/2);
        xmaxidx = xmaxidx + floor(dd/2);
    end

    speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);

    pl = imagesc(log10(speed_shelf));
    set(pl, 'AlphaData', ~isnan(speed_shelf));
    hold on
    inshelf = shelf.in;

    % add front and grounding line points
    isfrontshelf = data.isfront;
    isfrontshelf(~(shelf.in)) = 0;
    isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
    a = isfrontshelf;
    [x,y] = find(a);
    scatter(y,x,6,front_clr, 'filled')
    set(gca,'YDir','rev')

    isglshelf = data.isgl;
    isglshelf(~(shelf.in)) = 0;
    isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
    a = isglshelf;
    [x,y] = find(a);

    scatter(y,x,6,gl_clr, 'filled')
    set(gca,'YDir','rev')

    ax(i).XTick = [];
    ax(i).YTick = [];
    axis equal
    ax(i).Visible = 'off';
    %c(i) = colorbar(ax(i));
    clim([0,4])
    fig(i).Position(3:4) = [560, 420];
    colormap(cmap)

    % add a 50km line
    plot(0:100, ones(1,101), 'k', 'linewidth', 2); %500m resolution

    fprintf(shelves(i))
    fprintf(': number of grounding line points is %.0f \n', sum(sum(isglshelf)));


end

%% Pine Island
fig = figure(1); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/PIG.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 20; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
%p = spy(isfrontshelf);


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
spy(isglshelf, 'k');


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, zeros(1,101), 'k', 'linewidth', 2)

%add the dimensions of the box
box_y = [217, 268];
box_x = [207, 253];

if ~zoom
    plot([box_x(1), box_x(2)], [box_y(1), box_y(1)], 'k--', 'linewidth', 1.2);
    plot([box_x(1), box_x(2)], [box_y(2), box_y(2)], 'k--', 'linewidth', 1.2);
    plot([box_x(1), box_x(1)], [box_y(1), box_y(2)], 'k--', 'linewidth', 1.2);
    plot([box_x(2), box_x(2)], [box_y(1), box_y(2)], 'k--', 'linewidth', 1.2);
end

fprintf('number of PIG grounding line points is %.0f \n', sum(sum(isglshelf)));

if zoom % do the zoomed in thing

    uu = u_shelf(logical(isglshelf));
    vv = v_shelf(logical(isglshelf));
    [x,y] = find(isglshelf);
    quiver(y,x,vv,uu, 'LineWidth',1.25, 'color', 0.75*[1,1,1]);

    % put 2km on
    plot([box_x(1), box_x(1)+4],[box_y(2), box_y(2)], 'k', 'linewidth', 2)
    spy(isglshelf, 'k',15 );
    ax(1).XLim = box_x;
    ax(1).YLim = box_y;

    fprintf('maximum velocity is %.3f m/yr \n', sqrt(max(uu.^2 + vv.^2)));

end

%% Ross
fig = figure(2); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/Ross.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 50; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
%p = spy(isfrontshelf);


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
spy(isglshelf, 'k');


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, zeros(1,101), 'k', 'linewidth', 2); %500m resolution

fprintf('number of Ross grounding line points is %.0f \n', sum(sum(isglshelf)));

%% Ronne
fig = figure(3); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/Ronne.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 50; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
%p = spy(isfrontshelf);


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
spy(isglshelf, 'k');


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, zeros(1,101), 'k', 'linewidth', 2); %500m resolution

fprintf('number of Ronne grounding line points is %.0f \n', sum(sum(isglshelf)));

%% Amery
fig = figure(4); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/Amery.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 50; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
%p = spy(isfrontshelf);


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
spy(isglshelf, 'k');


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, zeros(1,101), 'k', 'linewidth', 2); %500m resolution

fprintf('number of Amery grounding line points is %.0f \n', sum(sum(isglshelf)));


%% Filchner
fig = figure(5); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/Filchner.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 50; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
%p = spy(isfrontshelf);


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
spy(isglshelf, 'k');


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, zeros(1,101), 'k', 'linewidth', 2); %500m resolution

fprintf('number of Filchner grounding line points is %.0f \n', sum(sum(isglshelf)));

%% Larsen
fig = figure(6); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/Larsen.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 50; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
%p = spy(isfrontshelf);


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
spy(isglshelf, 'k');


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, zeros(1,101), 'k', 'linewidth', 2); %500m resolution

fprintf('number of Larsen grounding line points is %.0f \n', sum(sum(isglshelf)));

%% Getz
fig = figure(7); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/Getz.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 50; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
%p = spy(isfrontshelf);


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
spy(isglshelf, 'k');


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, zeros(1,101), 'k', 'linewidth', 2); %500m resolution

fprintf('number of Getz grounding line points is %.0f \n', sum(sum(isglshelf)));

%% PSK
fig = figure(8); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/PopeSmithKohler.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 50; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
a = isfrontshelf;
[x,y] = find(a);
scatter(y,x,6,front_clr, 'filled')
set(gca,'YDir','rev')


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
a = isglshelf;
[x,y] = find(a);

scatter(y,x,6,gl_clr, 'filled')
set(gca,'YDir','rev')


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, ones(1,101), 'k', 'linewidth', 2); %500m resolution

fprintf('number of PSK grounding line points is %.0f \n', sum(sum(isglshelf)));

%% Thwaites
fig = figure(7); clf;
ax(1) = gca;

fname = strcat("data-for-figures/shelves/Thwaites.mat");
shelf = load(fname);

[rId, cId] = find(shelf.in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
np = 50; % number of padding points
xminidx = xminidx - np;
xmaxindx = xmaxidx + np;
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
yminidx = yminidx - np;
ymaxindx = ymaxidx + np;


dy = ymaxidx - yminidx;
dx = xmaxidx - xminidx;
dd = abs(dx - dy);
if dy < dx %adjust so square
    yminidx = yminidx - floor(dd/2);
    ymaxidx = ymaxidx + floor(dd/2);
else
    xminidx = xminidx - floor(dd/2);
    xmaxidx = xmaxidx + floor(dd/2);
end

speed_shelf = speed(xminidx:xmaxidx, yminidx:ymaxidx);
u_shelf = u(xminidx:xmaxidx, yminidx:ymaxidx);
v_shelf = v(xminidx:xmaxidx, yminidx:ymaxidx);

pl = imagesc(log10(speed_shelf));
set(pl, 'AlphaData', ~isnan(speed_shelf));
hold on
inshelf = shelf.in;
inshelf = inshelf(xminidx:xmaxidx, yminidx:ymaxidx);

%contour(inshelf, [0.5, 0.5], 'k', 'linewidth', 1.5) %add the grounding line and ice front

% add front and grounding line points
isfrontshelf = data.isfront;
isfrontshelf(~(shelf.in)) = 0;
isfrontshelf = isfrontshelf(xminidx:xmaxidx, yminidx:ymaxidx);
a = isfrontshelf;
[x,y] = find(a);
clr = a(a~=0);
clr = 0.75*[1,1,1];
scatter(y,x,6,clr, 'filled')
set(gca,'YDir','rev')


isglshelf = data.isgl;
isglshelf(~(shelf.in)) = 0;
isglshelf = isglshelf(xminidx:xmaxidx, yminidx:ymaxidx);
a = isglshelf;
[x,y] = find(a);
clr = a(a~=0);
clr = [0,0,0];
scatter(y,x,6,clr, 'filled')
set(gca,'YDir','rev')


ax(1).XTick = [];
ax(1).YTick = [];
axis equal
ax(1).Visible = 'off';
%c(i) = colorbar(ax(i));
clim([0,4])
fig(1).Position(3:4) = [560, 420];
colormap(cmap)

% add a 50km line
plot(0:100, ones(1,101), 'k', 'linewidth', 2); %500m resolution

fprintf('number of Thwaites grounding line points is %.0f \n', sum(sum(isglshelf)));


%% Make the colorbar
figure(10); clf;
c = colorbar;
colormap(cmap);
c.Ticks = 0:0.25:1;
c.TickLabels = {'10^0', '10^1', '10^2', '10^3', '10^4'};
axs = gca; axs.Visible = 'off';