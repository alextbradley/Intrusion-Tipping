% Generate files for the ice shelves

% Generate indices and lengthscales for some ice shelves

folder = 'shelves/';
fname = 'PIGfast'; %name of the ice shelf
fname = strcat(folder, fname);
%f = load('Antarctica-data.mat');
%%

figure(1); clf;
spy(f.isfront, 'b')
hold on 
spy(f.isgl, 'r')
%v = sqrt((f.vx).^2 + (f.vy).^2);
%h = imagesc(f.x, f.y, f.mask);
set(gca, 'YDir', 'normal')
axis equal
title('zoom and click around ice shelf (blue)')
zoom on;
pause
%waitfor(gcf, 'CurrentCharacter', char(13)) %what until we hit return
zoom reset
zoom off; % to escape the zoom mode

[xf,yf] = ginput();
hold on; plot(xf,yf,'b');

% title('click around the ice shelf grounding line (red)')
% [xg,yg] = ginput();
% plot(xg,yg, 'ro-')

%%
% [xx,yy] = meshgrid(f.y, f.x);
% in_shelf = inpolygon(xx, yy, xf,yf);

xidx = 1:10229;
yidx = 1:10941;
[xx,yy] = meshgrid(xidx,yidx);
in_shelf = inpolygon(xx, yy, xf,yf);

isgl_shelf = f.isgl & in_shelf;
isfront_shelf = f.isfront & in_shelf;

isfront = (f.isfront);
isfront(~in_shelf) = 0;
isgl = f.isgl;
isgl(~in_shelf) = 0;
figure(2); clf; spy(isfront);
hold on
spy(isgl, 'r')
set(gca, 'YDir', 'normal')

%% save 
bedmap_shelf = (f.mask == 3);
in = in_shelf & bedmap_shelf; 
fname = strcat(fname, '.mat');
save(fname, 'in ')

