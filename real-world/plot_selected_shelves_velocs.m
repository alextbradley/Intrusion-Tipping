% Make a plot of the velocities of selected shelves. Add grounding line
% points in red and ice front points in black. Add ocean forcing around the
% edges.

% Filenmaes
fnames = ["Amery", "Filchner", "Getz", "Larsen","PIGfast", "PopeSmithKohler", "Ronne", "Ross", "Thwaites"];

% Get the shelf data
f = load('Antarctica-data.mat');
mask = f.mask;
bedmap_shelf = mask == 3; %bedmap shelf points
v = sqrt(f.vx.^2 + f.vy.^2); %ice velocities

% Loop over file names, add shelf indices to each
shelves = zeros(size(v));
for i = 1:length(fnames)
    fname = strcat('shelves/', fnames(i), '.mat');
    in = load(fname);
    in = in.in; %target indices
    %in = in & bedmap_shelf; 
    shelves = shelves | in; %
    %save(fname, "in");
end

%remove non shelf points
v(shelves == 0) = nan;


%% make plot
figure(1); clf; 
hh = imagesc(f.tf_max);
colormap(gca, cmocean('speed'));
%c2 = colorbar;
clim([0, 4]);
shg
ax = gca;
axis equal
ax2 = axes();
ax2.Position = ax.Position;
h = imagesc(ax2, v);
shg
set(h, 'AlphaData', ~isnan(v)); %make nans see through
ax2.Visible = 'off';
colormap(ax2, cmocean('dense')); clim([0, 4000])
set(hh, 'AlphaData', ~isnan(f.tf_max));
axis equal
ax2.XLim = ax.XLim;
ax2.YLim = ax.YLim;

ax3 = axes;
ax3.Position = ax2.Position;
hold on
spy(isgl, 'r')
spy(isfront, 'k');
ax3.Visible = 'off';
axis equal
ax3.XLim = ax.XLim;
ax3.YLim = ax.YLim;

%% 