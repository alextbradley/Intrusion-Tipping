% Make supplementary figure 10 of the manuscript, showing the intrusion
% length as a function of time for different values of the tidal velocity.
%
% 06/11/23, ATB (aleey@bas.ac.uk). MIT licence
%
%% Preliminaries
addpath('functions');

% figure setup
basecolor = [0, 0,1; %cold
    1,0, 0]; %warm

colmap = nan(2,5,3);
cols = linspecer(5);
for i = 1:2
    for j = 1:5
        %colmap(i,j,:) = basecolor(i,:)*(1-j/10);
        colmap(i,j,:) = cols(j,:);

    end
end


% load the data
data = load('data-for-figures/supfigure10_data.mat', 'sol');
sol = data.sol;

%parameters
lengthscale = 25; %x lengthscale

%% Make the plot
figure(1); clf; hold on; box on;

for i = 1:2
    for j = 1:5

        intrusion_length = nan(1,length(sol(i,j).t));
        t = nan(1,length(sol(i,j).t));


        for iT = 1:(length(sol(i,j).t)-1)
            %if we're in second case, add first row in low alpha

            idx = find(sol(i,j).h(iT,:)-sol(i,j).h1(iT,:) == 0, 1, 'First');
            intrusion_length(iT) = -sol(i,j).x(idx)*lengthscale;
            t(iT) = sol(i,j).t(iT);
        end
        if i==2
            plot(t, intrusion_length, 'linewidth', 1.75, 'Color',colmap(i,j,:));
            legendinfo{j} = ['$u_t/u_{\inf} = ' num2str(sol(i,j).pp.u_tidal) '$'];
        else
            plot(t, intrusion_length, '--', 'linewidth', 1.75, 'Color',colmap(i,j,:), 'HandleVisibility','off');

        end


        %       
%   drawnow
%         pause

    end
end

% add the figure 2 results
fig2data = load('data-for-figures/figure2_data.mat');
fig2sol = fig2data.sol;


for i = 1:2

        intrusion_length = nan(1,length(fig2sol(i).t));
        t = nan(1,length(fig2sol(i).t));


        for iT = 1:(length(fig2sol(i).t)-1)
            %if we're in second case, add first row in low alpha

            idx = find(fig2sol(i).h(iT,:)-fig2sol(i).h1(iT,:) == 0, 1, 'First');
            intrusion_length(iT) = -fig2sol(i).x(idx)*lengthscale;
            t(iT) = fig2sol(i).t(iT);
        end
        if i==2
            plot(t, intrusion_length, 'linewidth', 1.5, 'Color','k');
        else
            plot(t, intrusion_length, '--', 'linewidth', 1.5, 'Color','k');

        end
        %       
%   drawnow
%         pause

end


xlabel('time (days)');
ylabel('intrusion distance (m)');
set(gca, 'YScale', 'log');
xlim([0,5])
legend(legendinfo, 'Interpreter','latex', 'location', 'northwest');
ax = gca;
ax.FontSize = 14;
