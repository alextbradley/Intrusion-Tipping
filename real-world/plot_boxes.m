% make boxes for the shelves
folder  = 'shelves/';
shelves = ["Amery", "Filchner", "Getz", "Larsen", "PIGfast", "PopeSmithKohler", "Ronne", "Ross", "Thwaites"];

sz = size(shelves);
mean_tf_max = nan(sz);
std_tf_max  = nan(sz);
mean_slope  = nan(sz);
std_slope  = nan(sz);
mean_velocs = nan(sz);
std_velocs  = nan(sz);

for i = 1:length(shelves)
    fname = strcat(folder, shelves(i), '.mat');
    f = load(fname);

    %thermal forcing
    tf_max = f.tf_shelf;
    tf_max = tf_max(~isnan(tf_max));
    mean_tf_max(i) = mean(tf_max);
    std_tf_max(i) = std(tf_max);

    %slope
    slope = f.slope;
    slope = slope(~isnan(slope));
    mean_slope(i)= mean(slope);
    std_slope(i) = std(slope);

    %velocity
    velocs = f.velocs;
    velocs = velocs(~isnan(velocs));
    mean_velocs(i)= mean(velocs);
    std_velocs(i) = std(velocs);

end

%% Compute dimensionless quantities
L = 335000; 
c = 3974;
St = 5.9e-4;
Cd = 1e-2;
uinf = 0.01; % 1 cm/s
secs_per_year = 365*24*60^2; %ice velocities are in m/a

dT = mean_tf_max / (L/c)  * St / Cd * uinf ./ mean_velocs * secs_per_year;
dT_std_pos =  (mean_tf_max + 0.5*std_tf_max) / (L/c)  * St / Cd * uinf ./ (mean_velocs - 0.5*std_tf_max) * secs_per_year; %one std higher
dT_std_neg =  (mean_tf_max - 0.5*std_tf_max) / (L/c)  * St / Cd * uinf ./ (mean_velocs + 0.5*std_tf_max) * secs_per_year; %one std lower

S = tan(mean_slope)/Cd;
S_std_pos = tan(mean_slope + 0.5*std_slope)/Cd;
S_std_neg = tan(mean_slope - 0.5*std_slope)/Cd;


% plot as boxes
%figure(1);clf; hold on
colmap = lines(max(sz));
for i = 1:max(sz)
    xf = [S_std_neg(i), S_std_pos(i), S_std_pos(i), S_std_neg(i)];
    yf = [dT_std_neg(i), dT_std_neg(i), dT_std_pos(i), dT_std_pos(i)];
   % f = fill(xf, yf, colmap(i,:), 'FaceAlpha', 0.2, 'LineWidth', 1, 'EdgeColor', colmap(i,:));
    plot(S(i), dT(i), 'o', 'markerfacecolor', colmap(i,:), 'markeredgecolor', 'k')
    text(S(i)+ 0.1, dT(i), shelves(i))
end

xlabel('S');
ylabel('\Delta T')
%xlim([-0.4, 0.4])
%ylim([0,10])
