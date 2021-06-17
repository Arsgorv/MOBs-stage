function [Binsizes] = Theta_binsize(nmice, experiment, def_binsize)
% 
% This function calculate dynamical binsize based on the Theta-rythm/Delta-rythm ratio.
% INPUT
%       nmice               array of mice. Example: nmice = [905 906 911 994 1161 1162 1168 1182 1186];
%       experiment          type of the experiment. Currently only 'PAG' is available.
% OUTPUT
%       Binsizes
%       ThetaEpoch_nmice
%       S
%       AllPeaks
%
% By Arsenii Goriachenkov, MOBS team, Paris,
% 30/05/2021
% github.com/arsgorv

% Перевод tsd в сек: 1 сек = 10 000 tsd. LFP - в tsd Время в секундах порядка 2е4, в тсд 2е8, на выходе findpeaks пока что в 2е3

%% Temporary parameters
AddMyPaths_Arsenii;
experiment = 'PAG';
def_binsize = 1000; % 100 ms = 0.1 sec = 1000 ts
nmice = [905 906 911 994 1161 1162 1168 1182 1186];

%% Parameters
if strcmp(experiment, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(experiment, 'MFB')
    fetchpaths = 'StimMFBWake';
end

Options.Fs = 1250; % ???? sampling rate of LFP
Options.TimeLim = 0; % ???? in seconds, minimal distance between two minimums or maximums. Used to be 0.08
Options.FilBand = [5 12];
Options.std = [0.4 0.1]; % ???? std limits for first and second round of peak

%% UNUSED: Parameters to calculate sparsed LFP
% Options.WVDownsample = 10; % ???? Зачем мы разрежаем LFP и идем с шагом 10?
% Options_1 = Options;
% Options_1.Fs = Options.Fs/Options.WVDownsample;
% tps = Range((LFP));
% vals = Data((LFP));
% LFPdowns = tsd(tps(1:Options.WVDownsample:end), vals(1:Options.WVDownsample:end));
% AllLFP.H = LFPdowns;

%% Allocate memory
ch_hpc = cell(1, length(nmice));
LFP = cell(1, length(nmice));
AllPeaks = cell(1, length(nmice));
Binsizes = cell(1, length(nmice));
Hz_Binsizes = cell(1, length(nmice));
S = cell(1, length(nmice));

%% Load the data

% Get paths of each individual mouse
Dir = PathForExperimentsERC_Arsenii(fetchpaths);
Dir = RestrictPathForExperiment(Dir, 'nMice', nmice);

for imouse = 1:length(Dir.path)
    try
        ch_hpc{imouse} = load([Dir.path{imouse}{1} 'ChannelsToAnalyse\dHPC_deep.mat']); % Number of the channel used
    catch
        ch_hpc{imouse} = load([Dir.path{imouse}{1} 'ChannelsToAnalyse\dHPC_rip.mat']);
    end
    LFP{imouse} = load([Dir.path{imouse}{1} 'LFPData\LFP' num2str(ch_hpc{1, imouse}.channel) '.mat']); % LFP from a channel ch_hpc
%     figure, plot(Range(LFP{1, imouse}.LFP, 'sec'), Data(LFP{1, imouse}.LFP))
end

%% Calculate spectrum

for imouse = 1:length(Dir.path)
    LowSpectrum_AG([strip(Dir.path{imouse}{1}, '\') filesep], ch_hpc{1, imouse}.channel,'H', imouse, Dir); % You can change spectrogramm windown size in 19th row (69 row)
    S{imouse} = load([Dir.path{imouse}{1} 'H_Low_Spectrum.mat']); % (1) - time x freq (2) - time (now window is - 200 ms). (3) - freq in Hz.
end

%% Define theta intervalset

for imouse = 1:length(Dir.path)

    % Cut frequencies to Theta and Delta range
    id_start_T = find(S{1, imouse}.Spectro{3} >= 6, 1, 'first');
    id_end_T = find(S{1, imouse}.Spectro{3} <= 10, 1, 'last');
    
    id_start_D = find(S{1, imouse}.Spectro{3} >= 2, 1, 'first');
    id_end_D = find(S{1, imouse}.Spectro{3} <= 5, 1, 'last');
    
    Theta_D = S{1, imouse}.Spectro{1}(:,id_start_T:id_end_T);
    Delta_D = S{1, imouse}.Spectro{1}(:,id_start_D:id_end_D);
    
    Theta = mean(Theta_D,2);
    Delta = mean(Delta_D,2);
    Ratio = Theta./Delta;
    
    % Choose the threshold to define theta intervalset
    figure('Units','normalized','Position',[0 0 1 1]);
    subplot(211);
    
    imagesc(S{1, imouse}.Spectro{2}, S{1, imouse}.Spectro{3}, S{1, imouse}.Spectro{1}');
    axis xy;
    caxis([0 8e4]);
    title(['Spectro of mouse #' num2str(nmice(imouse))]);
    
    i = 0;
    while i == 0
        % Enter the caxis upper value
        prompt = {'Enter the caxis upper value'};
        definput = {'8e4'};
        dlgtitle = 'Input';
        dims = [1 35];
        caxis_v = inputdlg(prompt, dlgtitle, dims, definput);
        caxis_v = str2double(caxis_v{1});
        caxis([0 caxis_v]);
        
        % Are you satisfied?
        prompt = {'Are you satisfied? (y/n)'};
        sat = inputdlg(prompt);
        if sat{1} == 'n'
            i = 0;
        else
            i = 1;
        end
        clear sat prompt
    end
    
    subplot(212);
    plot(S{1, imouse}.Spectro{2},runmean(Ratio, 60));
    xlim([0 S{1, imouse}.Spectro{2}(end)]);
    
    i = 0;
    while i == 0
        [~, theta_threshold] = ginput(1);
        y_line = line(xlim, [theta_threshold theta_threshold], 'Color', 'r');
        
        % Are you satisfied with the threshold?
        prompt = {'Are you satisfied? (y/n)'};
        sat = inputdlg(prompt);
        if sat{1} == 'n'
            i = 0;
        else
            i = 1;
            close
        end
        delete(y_line);
    end

    Ratio_tsd = tsd((S{1, imouse}.Spectro{2}*1e4)', runmean(Ratio, 60));
    IntervalTheta = thresholdIntervals(Ratio_tsd, theta_threshold, 'Direction', 'Above');
    merge_value = 2; % in seconds. 1s = 10 000 tsd
    drop_value = 3; % in seconds. 1s = 10 000 tsd
    ThetaEpoch = mergeCloseIntervals(IntervalTheta, merge_value*1e4);
    ThetaEpoch = dropShortIntervals(ThetaEpoch, drop_value*1e4);
    ThetaEpoch_nmice{imouse} = ThetaEpoch;
    
    % Save the result
    save(['E:\ERC_data\M' num2str(nmice(imouse)) '\ThetaEpoch.mat'], 'ThetaEpoch', 'merge_value', 'drop_value', 'theta_threshold');
    save(['C:\Users\Arsenii Goriachenkov\Documents\MATLAB\Theta\ThetaEpoch_nmice.mat'], 'ThetaEpoch_nmice');
end

%% Sort theta and non-theta epochs and create the sorted list of intervals

% Create Theta and Non-Theta intervalsets
for imouse = 1:length(Dir.path)  
    global_time{imouse} = Range(LFP{imouse}.LFP); %в ts (10000 ts = 1sec)
    GlobalEpoch{imouse} = intervalSet(global_time{imouse}(1), global_time{imouse}(end));
    non_ThetaEpoch_nmice{imouse} = GlobalEpoch{imouse} - ThetaEpoch_nmice{imouse};
end

% Create arrays of Theta and Non-Theta Start-End intervals
for imouse = 1:length(Dir.path)
    Start_theta{imouse} = Start(ThetaEpoch_nmice{imouse});
    End_theta{imouse} = End(ThetaEpoch_nmice{imouse});
    Theta_array{imouse} = [Start_theta{imouse} End_theta{imouse}];

    Start_non_theta{imouse} = Start(non_ThetaEpoch_nmice{imouse});
    End_non_theta{imouse} = End(non_ThetaEpoch_nmice{imouse});
    non_Theta_array{imouse} = [Start_non_theta{imouse} End_non_theta{imouse}];
end

% Merge and sort
for imouse = 1:length(Dir.path)
    
    i = 1;
    j = 1;
    Map_TimeLine{imouse} = nan(length(Start_non_theta{imouse})+length(Start_theta{imouse}), 3); % in ts (10000 ts = 1 sec)
    
    if Start_theta{imouse}(1) < Start_non_theta{imouse}(1) % if Theta - first to come
        for i = 1:length(Start_non_theta{imouse})
            Map_TimeLine{imouse}(j, 1) = Start_theta{imouse}(i);
            Map_TimeLine{imouse}(j, 2) = End_theta{imouse}(i);
            Map_TimeLine{imouse}(j, 3) = 1; % means - theta interval
            
            Map_TimeLine{imouse}(j+1, 1) = Start_non_theta{imouse}(i);
            Map_TimeLine{imouse}(j+1, 2) = End_non_theta{imouse}(i);
            Map_TimeLine{imouse}(j+1, 3) = -1; % means - non theta interval
            
            j = j + 2;
        end
        Map_TimeLine{imouse}(j, 1) = Start_theta{imouse}(i+1);
        Map_TimeLine{imouse}(j, 2) = End_theta{imouse}(i+1);
        Map_TimeLine{imouse}(j, 3) = 1; % means - theta interval
        
    else % If non-Theta - first to come
        for i = 1:length(Start_theta{imouse})
            Map_TimeLine{imouse}(j, 1) = Start_non_theta{imouse}(i);
            Map_TimeLine{imouse}(j, 2) = End_non_theta{imouse}(i);
            Map_TimeLine{imouse}(j, 3) = -1; % means - non theta interval
             
            Map_TimeLine{imouse}(j+1, 1) = Start_theta{imouse}(i);
            Map_TimeLine{imouse}(j+1, 2) = End_theta{imouse}(i);
            Map_TimeLine{imouse}(j+1, 3) = 1; % means - theta interval
            
            j = j + 2;
        end
        Map_TimeLine{imouse}(j, 1) = Start_non_theta{imouse}(i+1);
        Map_TimeLine{imouse}(j, 2) = End_non_theta{imouse}(i+1);
        Map_TimeLine{imouse}(j, 3) = -1; % means - non theta interval
    end
    save(['E:\ERC_data\M' num2str(nmice(imouse)) '\Map_TimeLine.mat'], 'Map_TimeLine');
end


%% Caculate peaks to define phase of theta-rythm

% Calculate times at which curve has maximum and minimum peaks
for imouse = 1:length(Dir.path)
    AllPeaks = FindPeaksForFrequency_AG(LFP{1, imouse}.LFP, Options, 0);
    AllPeaks_nmice{imouse} = AllPeaks;
    save(['E:\ERC_data\M' num2str(nmice(imouse)) '\AllPeaks.mat'], 'AllPeaks');
end

% Filter AllPeaks by maximums and restrict peaks to theta-epoch
for imouse = 1:length(Dir.path)
    filmin_AllPeaks_nmice{imouse} = ts(AllPeaks_nmice{imouse}(find(AllPeaks_nmice{imouse}(:, 2) == 1))*1e4); % filtered by maxs
    Theta_Peaks{imouse} = Restrict(filmin_AllPeaks_nmice{imouse}, ThetaEpoch_nmice{imouse});
end

save(['C:\Users\Arsenii Goriachenkov\Documents\MATLAB\Theta\Theta_Peaks.mat'], 'Theta_Peaks');
save(['C:\Users\Arsenii Goriachenkov\Documents\MATLAB\Theta\Map_TimeLine.mat'], 'Map_TimeLine');

%% Trash
%
% [maxpower, id_mf] = max(Theta_D, [], 2); % на каждом временном окне нашли значение частоты с максимальной мощностью
% Theta_freq = S{3}(id_start_D:id_end_D);
% maxfreq = arrayfun(@(x) Theta_freq(id_mf(x)), 1:length(id_mf)); % наиболее представленная частота
% figure, plot(maxfreq) %визуализация динамики частоты
% %% Фильтрация LFP
% FilLFP_hpc = FilterLFP(LFP,[6 13],1024);
% figure, plot(Range(FilLFP_hpc, 's'), Data(FilLFP_hpc))

end