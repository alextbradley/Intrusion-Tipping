% Put the following data onto the 1km grid defined by Adusumilli et al:

% Ice thickness (used to determine grounding line position, from Bedmachine
% V3)
% Thermal forcing (downsampled from 5km data of Adusumilli et al 2020)
% Ice velocity (from ITS_LIVE)
% Bed depth (from Bedmachine V3)

%% define the grid
x = h5read('bb0448974g_3_1.h5', '/x');
y = h5read('bb0448974g_3_1.h5', '/y'); 
melt = h5read('bb0448974g_3_1.h5', '/w_b'); %load the melt for comparison
[xx,yy] = meshgrid(x,y);
xx = xx';
yy = yy';

%% Ice thickness, bed depth and mask from bedmachine
H = interpBedmachineAntarctica(xx,yy,'thickness');
B = interpBedmachineAntarctica(xx,yy,'bed');
mask = interpBedmachineAntarctica(xx,yy,'mask');

%% Thermal forcing 
% get the raw, 5km resolution data
tf_max_5km = h5read('max_thermal_forcing_200-800.h5', '/tf_max');
x_tf = h5read('max_thermal_forcing_200-800.h5', '/x');
y_tf = h5read('max_thermal_forcing_200-800.h5', '/y');

% make a scatteredinterpolant
[xx_tf,yy_tf] = meshgrid(x_tf,y_tf);
xx_tf = xx_tf';
yy_tf = yy_tf';
F = scatteredInterpolant(xx_tf(:), yy_tf(:), tf_max_5km(:));

% evaluate scattered interpolant on our grid
tf_max = F(xx, yy);

%% compute masks of ice front and grounding line 
sz = size(H);
isedge = zeros(sz);
isgl = zeros(sz);
isfront = zeros(sz);
tic
for i = 2:sz(1)-1
    for j = 2:sz(2)-1

        if (mask(i,j) == 3 && (any(mask(i-1:i+1,j) ~= 3) || any(mask(i,j-1:j+1) ~= 3) ))   %if we're in the shelf and a point next to us isn't
            isedge(i,j) = 1;

            %determine whether its a gl point or not
            if (all(H(i-1:i+1,j)) && all(H(i,j-1:j+1))) %&& H(i,j) > 0.95*abs(B(i,j))*1028/918) %if we have thickness at neighbouring pts
                isgl(i,j) = 1;
            elseif H(i,j) < 0.95*abs(B(i,j))*1028/918 %add a thickness based threshold
                isfront(i,j) = 1;
            end
        end
    end
end
toc

%% get the velocities from its_live
addpath('./ITS_LIVE_tools/Antartic-Mapping-Tools/')
addpath('./ITS_LIVE_tools/ITS_LIVE-tools/')
variable = 'vx'; vx = itslive_interp(variable,xx,yy);
variable = 'vy'; vy = itslive_interp(variable,xx,yy);
v = sqrt(vx.^2 + vy.^2); 


%% save the data
%save('Antarctica-data.mat', 'x','y', 'tf_max', 'H', 'B', 'mask', 'isedge', 'isgl', 'isfront')
