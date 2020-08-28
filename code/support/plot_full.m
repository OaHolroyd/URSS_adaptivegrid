%% PLOT FULL
% Plots a single output
% INPUTS:  Tout - output times
%          Xout - non-uniform grids at Tout
%          Uout - interface height at Xout/Tout
%          Fout - control/forcing term at Xout/Tout
%          Wout - weight function

sf = 2^2; % number of steps to jump between frames
p_time = 0.0; % pause between frames

start_index = 1; % start plotting here
end_index = length(Tout); % stop plotting here

wsf = 1; % scale up W by wsf if it's too small to see

%% Initial Plot

% create figure
fig = figure(1);

hold on
uplot = plot(Xout{start_index},Uout{start_index}-1,'blue');
fplot = plot(Xout{start_index},Fout{start_index},'red');
wplot = plot(0,0,'cyan');
pplot = plot(0,0,'magenta');
hold off

lgd = legend({['non-uniform (N = ',num2str(length(Xout{start_index})),')'],'forcing term','weight function','padded weight'});
xlabel('$x$')
ylabel('interfacial height, forcing term, grid weighting')

% set axis limits, scale and position
ylim([-0.3,0.7])
xlim([0,L])
pbaspect([6,3,1])

ax = gca;
ax.XTick = [Xout{1};L];
ax.XTickLabel = {};
ax.XGrid = 'on';

% this probably needs changing depending on if/where your second screen is
fig.Position = [1441, -179, 1920, 984];

pause(p_time);


%% Animated Plot

for i = (1+start_index):sf:end_index
    if sum(abs(Uout{i})>50000)>0
        disp(['break at i = ',num2str(i), ' (t = ',num2str(Tout(i)),')'])
        return
    end
    
    ax.XTick = [Xout{i};L];
    ax.XTickLabel = {};
    ax.XGrid = 'on';
    
    fig.Name = strcat("t = ",num2str(Tout(i)));
    lgd = legend({['non-uniform (N = ',num2str(length(Xout{i})),')'],'forcing term','weight function','padded weight'});
    
    uplot.XData = Xout{i};
    uplot.YData = Uout{i}-1;
    fplot.XData = Xout{i};
    fplot.YData = Fout{i};
    wplot.XData = Xout{i};
    wplot.YData = wsf*Wout{i}; % this bit is probably what's going wrong if it fails to plot properly
    pplot.XData = Xout{i};
    pplot.YData = wsf*Pout{i};
    drawnow;
    pause(p_time);
end


%% Final plot

fig.Name = strcat("t = ",num2str(Tout(end_index))," (end)");
uplot.XData = Xout{end_index};
uplot.YData = Uout{end_index}-1;
fplot.XData = Xout{end_index};
fplot.YData = Fout{end_index};
wplot.XData = Xout{end_index};
wplot.YData = wsf*Wout{end_index}; % this bit is probably what's going wrong if it fails to plot properly
drawnow;











