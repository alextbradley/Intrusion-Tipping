% Generate the file 'all-Ant-data.mat', which contains...

%% Preliminaries
clear
addpath('functions');
addpath('data');
addpath('data/ITS_LIVE_tools/Antartic-Mapping-Tools/')
addpath('data/ITS_LIVE_tools/ITS_LIVE-tools/')
%% define the grid
x = h5read('bb0448974g_3_1.h5', '/x');
y = h5read('bb0448974g_3_1.h5', '/y'); %Adusumilli et al. 2020 co-ordinates
[xx,yy] = meshgrid(x,y);
xx = xx';
yy = yy';

%% Get ice velocities from ITS_LIVE
variable = 'vx'; vx = itslive_interp(variable,xx,yy);
variable = 'vy'; vy = itslive_interp(variable,xx,yy);

%% Get ice thickness, bed and mask from bedmachine 
H = interpBedmachineAntarctica(xx,yy,'thickness');
B = interpBedmachineAntarctica(xx,yy,'bed');
mask = interpBedmachineAntarctica(xx,yy,'mask');

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


