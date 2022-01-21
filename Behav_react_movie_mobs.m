function vidfile = Behav_react_movie_mobs(experiment, nmice)
%
% The function creates a movie of a mouse moving along its trajectory in ERC experiment, while on other 
% subplots you may observe the change of the reactivation strength (PC score) of neuronal assemblies,
% which was calculated with the 'react_pca_ica.m' script in the same manner it was calculated in 
% Peyrache et al., 2010.
% 
% Note for MOBS team: 
% Mind that you have to pre-generate ['PCscore_' mouse_name '.mat'] file (which is replayGtsd output 
% of the 'react_pca_ica.m' script) and put it into the mouse folder.
%
% Written by Arsenii Goriachenkov, MOBS team, Paris
% Dec 2021
% github.com/arsgorv

%% Manage inputs
% Experiment. Currently only used for 'PAG'
if strcmp(experiment, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(experiment, 'MFB')
    fetchpaths = 'StimMFBWake';
elseif strcmp(experiment, 'Novel')
    fetchpaths = 'Novel';
end

%% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',nmice);
[numsessions, micenames] = CountNumSesionsERC(Dir);

clear fetchpaths

%% Pre-allocate data
b = cell(numsessions, 1);
r = cell(numsessions, 1);

%% Load data
cnt = 1;
for imouse = 1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        warning('');
        
        % Load behavioural data
        b{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],'SessionEpoch', 'CleanVtsd',...
            'CleanAlignedXtsd', 'CleanAlignedYtsd', 'TTLInfo');
        
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            b{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],'SessionEpoch', 'Vtsd',...
                'AlignedXtsd', 'AlignedYtsd', 'ZoneEpoch', 'TTLInfo');
        end
        
        % Load ripples data
        try
            r{cnt} = load([Dir.path{imouse}{isession} 'SWR.mat'],'ripples');
        catch
            r{cnt} = load([Dir.path{imouse}{isession} 'Ripples.mat'],'ripples');
        end
        cnt = cnt + 1;
        
        % Load PC score data
        load([Dir.path{imouse}{isession} 'PCscore_' num2str(nmice) '.mat'], 'replayGtsd');
    end
end
clear cnt warnMsg imouse warnId warning

%% Choose parameters
list_epochs = {'TestPre1', 'TestPre2', 'TestPre3', 'TestPre4',...
    'Cond1', 'Cond2', 'Cond3', 'Cond4',...
    'TestPost1', 'TestPost2', 'TestPost3', 'TestPost4',};
[id_epoch,~] = listdlg('PromptString', {'Choose the epoch to create a movie'}, 'ListString', list_epochs, 'SelectionMode', 'single');
epoch_name = list_epochs{id_epoch};
epoch = eval(['b{isession}.SessionEpoch.' epoch_name]);

clear list_epochs id_epoch

%% Prepare behavioural data
try
    target_epoch_x = Restrict(b{isession}.CleanAlignedXtsd, epoch);
    target_epoch_y = Restrict(b{isession}.CleanAlignedYtsd, epoch);
    speed_epoch = Restrict(b{isession}.CleanVtsd, epoch);
catch
    target_epoch_x = Restrict(b{isession}.AlignedXtsd, epoch);
    target_epoch_y = Restrict(b{isession}.AlignedYtsd, epoch);
    speed_epoch = Restrict(b{isession}.Vtsd, epoch);
end

% Find shock times
shock_times = Range(Restrict(ts(Start(b{isession}.TTLInfo.StimEpoch)), epoch));

%% Prepare PCscore data

replayGtsd = (replayGtsd{1})';

% Pre-allocate data
PCscore = cell(length(replayGtsd), 1);
PCscore_time = cell(length(replayGtsd), 1);
PCscore_data = cell(length(replayGtsd), 1);

for itemp = 1:length(replayGtsd)
    PCscore{itemp} = Restrict(replayGtsd{itemp}, epoch);
    PCscore_time{itemp} = Range(PCscore{itemp});
    PCscore_data{itemp} = Data(PCscore{itemp});
    
    % Synchronize PCscore and behav times
        % Options below depend on the binsize of the spike train (we take 25 ms or 40 Hz)
        % And sampling-size of the behaviour (which is typically 14 Hz ~ 70 ms)
    
    if length(target_epoch_x) > length(PCscore{itemp})
        PCscore{itemp} = Restrict(PCscore{itemp}, target_epoch_x, 'align', 'closest');
    else %It doesn`t really change with itemp. So I leave the itemp=length(replayGtsd) value cause it doesn`t matter.
        target_epoch_x = Restrict(target_epoch_x, PCscore{itemp}, 'align', 'closest');
        target_epoch_y = Restrict(target_epoch_y, PCscore{itemp}, 'align', 'closest');
        speed_epoch = Restrict(speed_epoch, PCscore{itemp}, 'align', 'closest');
    end
end

%% Prepare ripples data
ripples = Restrict(ts(r{isession}.ripples(:,2)*1e4), epoch);

%% Calculate mean-among-templates PC score

% Pre-allocate data
count_temp = cell(length(replayGtsd), 1);
mean_PCscore = cell(length(PCscore_data{itemp}), 1);

for t = 1:length(PCscore_data{itemp})
    count = 1;
    while count <= length(replayGtsd)
        count_temp{count} = PCscore_data{count}(t);
        count = count + 1;
    end
    mean_PCscore{t} = mean(cell2mat(count_temp));
end
mean_PC = cell2mat(mean_PCscore);

clear count t count_temp mean_PCscore itemp

%% Other parameters
maze = [0 0; 0 1; 1 1; 1 0; 0.63 0; 0.63 0.75; 0.35 0.75; 0.35 0; 0 0];
shockZone = [0 0; 0 0.43; 0.35 0.43; 0.35 0; 0 0]; 
safeZone = [1 0; 0.63 0; 0.63 0.43; 1 0.43; 1 0]; 

%% Convert data

x_traj = Data(target_epoch_x);
y_traj = Data(target_epoch_y);
speed = Data(speed_epoch);
times = Range(target_epoch_x);

% Necessary assumption to align behahvioural sampling rate and ripples. 
ripples_aligned = Range(Restrict(ts(times), ripples, 'align', 'closest'));
shocks_aligned = Range(Restrict(ts(times), ts(shock_times), 'align', 'closest'));

%% Prepare video 
figh = figure;
figh.Position = [100 10 1200 750];
ishock = 1;

x_shock = zeros(length(x_traj), 1);
y_shock = zeros(length(x_traj), 1);
movieVector = struct('cdata', cell(1, length(x_traj)), 'colormap', cell(1, length(x_traj)));

% Create a VideoWriter object and set properties
vidfile = VideoWriter([micenames{isession} '_' epoch_name]);

vidfile.FrameRate = 22;

% Open the VideoWriter object, write the movie, and close the file
open(vidfile);

%% Render Scenario
for k = 1:length(x_traj) 
    
    % Clear the figure to start with a blank slate
    clf
    
    % Extract data at the current time step
    x_k = x_traj(k);
    y_k = y_traj(k);
    
    % Plot the current position
    fig1 = subplot(3, 7, 1);
    plot(x_k, y_k, 'o', 'LineWidth', 5, 'MarkerSize', 5, 'MarkerEdgeColor', [0.1 0.5 0.4])
    
    % Plot the entire trajectory
    hold on
    plot(x_traj, y_traj, 'k-.', 'LineWidth', 0.5);
    
    % Add plotting options
    grid on
    xlabel('x')
    ylabel('y')
    title([micenames{isession} ' ' epoch_name '. Time = ' num2str((times(k) - times(1))/1e4) ' seconds'])
    
    % Plot the maze
    plot(maze(:,1),maze(:,2),'k','LineWidth', 1.5)
    plot(shockZone(:,1), shockZone(:,2), 'color', [0.9 0.3 0.3], 'LineWidth', 1.5)
    plot(safeZone(:,1), safeZone(:,2), 'color', [0.1 0.7 0.5], 'LineWidth', 1.5)
    
    % Plot shocks
    if strcmp(epoch_name, 'Cond1') || strcmp(epoch_name, 'Cond2') || strcmp(epoch_name, 'Cond3') || strcmp(epoch_name, 'Cond4')
        if nonzeros((times(k) == shocks_aligned)) == 1
            plot(x_k, y_k, '*', 'LineWidth', 3, 'MarkerSize', 4, 'MarkerEdgeColor', [1 0 0])
            x_shock(ishock) = x_k;
            y_shock(ishock) = y_k;
        end
        if x_shock(1) ~= 0
            plot(nonzeros(x_shock(:)), nonzeros(y_shock(:)), '*', 'LineWidth', 3, 'MarkerSize', 4, 'MarkerEdgeColor', [1 0 0])
            ishock = ishock + 1;
            x_shock = nonzeros(unique(x_shock));
            y_shock = nonzeros(unique(y_shock));
        end
    end
    
    hold off

    % Define the time window
    fig_times = ((times(k)-2*1e4):(times(k)+2*1e4))';
    x_axis = linspace(fig_times(1), fig_times(end), 5);    
    t = (x_axis - times(1))/1e4;
    
    % Find ripples in the current window
    rips_onfig = find(ripples_aligned < fig_times(end));
    
    % Find shocks in the current window
    shocks_onfig = find(shocks_aligned < fig_times(end));
    
    % Plot the speed dynamics
    fig2 = subplot(3, 7, 2);
    h1 = plot(times, speed);
    title('Speed, m/s');
    ylim([0, max(speed)]);
    xlim([fig_times(1), fig_times(end)]);
    xticks(x_axis);
    xticklabels({t(1),t(2), t(3), t(4), t(5)});
    line([x_axis(3) x_axis(3)], ylim, 'Color', [0 0 0], 'LineWidth', 1);

    grid on;
    set(h1, 'Color', [0.5 0.1 0.5]);
        
    % Plot the mean PC score dynamics
    fig3 = subplot(3, 7, 3);
    h2 = plot(times, mean_PC);
    title('Mean PC score');
    ylim([min(mean_PC), max(mean_PC)]);
    xlim([fig_times(1), fig_times(end)]);
    xticks(x_axis);
    xticklabels({t(1),t(2), t(3), t(4), t(5)});
    line([x_axis(3) x_axis(3)], ylim, 'Color', [0 0 0], 'LineWidth', 1);

    grid on;
    set(h2, 'Color', [0.6 0.8 0.1]);
    
    % Plot ripples
    for i = 1:length(rips_onfig)
        line([ripples_aligned(rips_onfig(i)) ripples_aligned(rips_onfig(i))], ylim, 'Color', [0 0.2 0.5], 'LineWidth', 1.5);
    end
    clear i    
        
    % Plot shocks
    for i = 1:length(shocks_onfig)
        line([shocks_aligned(shocks_onfig(i)) shocks_aligned(shocks_onfig(i))], ylim, 'Color', [0.8 0 0], 'LineWidth', 1.5);
    end
    clear i
    
    % Plot PCs dynamics
    y_pos = 0.825;
    
    for jtemp = 1:length(replayGtsd)
        subplot(3, 7, 4);
        h3 = plot(times, PCscore_data{jtemp});
        set(h3, 'Color', [0.1 0.5 0.5]);
        ylim([min(PCscore_data{jtemp}), max(PCscore_data{jtemp})]);
        xlim([fig_times(1), fig_times(end)]);
        xticks(x_axis);
        xticklabels({t(1),t(2), t(3), t(4), t(5)});
        line([x_axis(3) x_axis(3)], ylim, 'Color', [0 0 0], 'LineWidth', 1);        
        grid on;
        title(jtemp);
        set(gca, 'Position', [0.6 y_pos 0.37 0.1]);
        y_pos = y_pos - 0.15;
        
        % Plot ripples
        for i = 1:length(rips_onfig)
            line([ripples_aligned(rips_onfig(i)) ripples_aligned(rips_onfig(i))], ylim, 'Color', [0 0.2 0.5], 'LineWidth', 1.5);
        end
        clear i
        
        % Plot shocks
        for i = 1:length(shocks_onfig)
            line([shocks_aligned(shocks_onfig(i)) shocks_aligned(shocks_onfig(i))], ylim, 'Color', [0.8 0 0], 'LineWidth', 1.5);
        end
        clear i
    end
            
    % Make figure pretty
    set(fig1, 'Position', [0.05 0.4 0.37 0.5]);
    set(fig2, 'Position', [0.05 0.1 0.37 0.05]);
    set(fig3, 'Position', [0.05 0.23 0.37 0.07]);
    
    % Take a Snapshot. Save the frame
    movieVector(k) = getframe(gcf);
    writeVideo(vidfile, movieVector(k));
end

close(vidfile);

end
