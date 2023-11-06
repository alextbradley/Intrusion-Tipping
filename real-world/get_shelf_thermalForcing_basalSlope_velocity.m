% Compute arrays of basal slope, thermal forcing and grounding line
% velocity for the ice shelves. Takes the input dataset for the entire
% Antarctcia, and a shape file ("IN") stored within the ice shelf mat
% file. Fname specifies the name of the shelf.
%
% 05/11/23 ATB (aleey@bas.ac.uk), MIT licence
%


folder = '../data-for-figures/shelves/';
fname = 'PopeSmithKohler'; %name of the ice shelf
fname = strcat(folder, fname, '.mat');
in = load(fname, 'in'); %indices of potential shelf points
in = in.in; 
%f = load('Antarctica-data.mat');
bedmap_shelf = (f.mask == 3);
in_shelf = in & bedmap_shelf; %just make sure we're only taking hself points, but this shouldn't change file 

% get the grounding line velocites
vx = f.vx;
vy = f.vy;
v = sqrt((f.vx).^2 + (f.vy).^2);

%ice shelf indices
isgl = f.isgl; 
isgl_shelf = isgl & in_shelf; %grounding line points in this shelf
isfront = f.isfront; 
isfront_shelf = isfront & in_shelf; %ice front points in the shelf

% shelf properties
vx_gl_shelf = vx(isgl_shelf);
vy_gl_shelf = vy(isgl_shelf);
v_gl_shelf  = sqrt(vx_gl_shelf.^2 + vy_gl_shelf.^2);

figure(1); 
subplot(2,2,1);
histogram(v_gl_shelf);
xlabel('grounding line velocity');
ylabel('count')

%% thermal forcing
% get the thermal forcing file
tf = f.tf_max;

% loop over the shelf front points
[rId, cId] = find(isfront_shelf) ; %row, column of front points
tf_shelf = nan(size(rId));
for ir = 1:length(rId)
    %take the mean over a small square
    square = tf(rId(ir)-3:rId(ir)+3, cId(ir)-3:cId(ir)+3);
    square = square(~isnan(square));
    tf_shelf(ir) = mean(square);
end

subplot(2,2,2)
histogram(tf_shelf);
xlabel('thermal forcing');
ylabel('count')
perc_coverage_thermal_forcing = sum(~isnan(tf_shelf))/ length(tf_shelf) * 100;
%% basal slope

% loop over the shelf front points
[rId, cId] = find(isgl_shelf) ; %row, column of front points
slope_sgl = nan(size(rId));
dx = 500;
dy = 500;
bed = f.B;
slope = zeros(size(rId));
velocs = zeros(size(rId));
for ir = 1:length(rId)
    %compute gradients in the two directions
    ix = rId(ir);
    iy = cId(ir);
    bx = (bed(ix+1,iy)-bed(ix-1,iy))/2/dx;
    by = (bed(ix,iy+1)-bed(ix,iy-1))/2/dy; %centered fd

    %get the flow direction from velocities
    vx_pt = vx(ix,iy);
    vy_pt = vy(ix,iy);
    velocs(ir) = sqrt(vx_pt^2 + vy_pt^2);
    direction =  [vx_pt;vy_pt];
    direction = direction/norm(direction);
    slope(ir) = [bx, by]*direction;
end

subplot(2,2,3)
histogram(slope, 1000);
xlabel('slope');
ylabel('count');
xlim([-0.05, 0.05])

subplot(2,2,4)
histogram(velocs);
xlabel('gl velocities (alt method)');
ylabel('count')

%% save data
save(fname, 'slope', '-append');
save(fname, 'velocs', '-append');
save(fname, 'tf_shelf', '-append');
