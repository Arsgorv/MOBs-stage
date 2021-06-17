% Clean_PWC
%
% This function performs pairwise correlations on all overlapping
% and non-overlapping place cells of each particular mouse and plots the
% results
%
% INPUT
%   
%   nmice               array with mice numbers that will harvested from ERC
%                       PathForExperiments. Each mouse should contain
%                       PlaceCells structure in its SpikeData.mat
%   SpeedThresh         threshold (in cm/s) to calculate place fields
%                       (default = 4) (optional)
%   BinSize             binsize (in tsd units) used to creat spike-time
%                       histograms (default = 0.1*1e4). In RepeatOriginalPaper
%                       mode equals 10 (1 ms) (optional)
%   OverlapFactor       defines how many pixels in place fields of two neurons
%                       should overlap to be considered overlapping
%                       (default = 25) (optional)
%   RepeatOriginalPaper if true, tt repeats the analysis of Wilson&McNaughton, 1994, Science
%                       (bins spike trains with 1ms bin and calculate average cross-correlation
%                       100 ms around). If false, calculates Pearson's 
%                       correlation coefficient on binned spiked traines
%                       (default = false) (optional)
%   RestrictSleep       if not empty, restrict sleep sessions to this time
%                       (in tsd units) (default = 20*60*1e4) (optional)
%   SaveFigure          if true, saves figures in dropbox (default = false)(optional)
%   PlotExamples        if true, plots cross-correlation histograms for
%                       every pair (optional)
% 
% OUTPUT
% 
%  None
% 
% EXAMPLE
% 
%   PairwiseCorrelationsAcrossMice_WM1994_DB([797 798 1117]);
%   PairwiseCorrelationsAcrossMice_WM1994_DB([797 798 1117], 'SaveFigure', true, 'BinSize', 0.05*1e4);
%
% By Arsenii Goriachenkov, MOBS team, Paris,
% 05/03/2021
%% Paths and functions
AddMyPaths_Arsenii

%% Parameters
% Behavioral states
states = {'PreSleep', 'Task', 'CondMov', 'CondFreeze', 'PostSleep', 'PostTests'};
% Colors for plots
cols = {[0.2 0.2 0.2], [0.8 0.8 0.8]};

%% mice which participate in analysis

% nmice = [905 906 911 994 1161 1162 1168 1182 1186]; %Mice in analysis
nmice = [905 906 911 994 1161 1162 1168 1182 1186]; %Mice in analysis
% nmice = [797 798 828 861 882 905 906 911 912 977 979 994 1117 1124 1161 1162 1168 1182 1186]; % All mice
% nmice = [797 798 861 882 977 1117]; % Не пройден PC threshold .
% 979 has no PAG

%% Optional parameters
% Speed threshold (cm/s). Epochs with lower speed are not considered.
speed_thresh = 4;
% Place cells threshold
PCnum_thresh = 2;
% Sleep time restriction. If not, put []
SleepTimeToRestrict = 20*60*1e4;
% Parameters of cross-correlograms. Measured in tsd units.
binsize = 0.1*1e4;
% Number of overlaped pixels to define a place cells` overlap
overlapFactor = 25;
% Do you want to correlate firing histograms (false) or short latency cross-correlations (true - very slow)?
do_as_in_paper = false;
% Do you want to save the figures? 
savefig = true;
% Do you want to plot each of single cross correlation of overlapping pairs?
plotexamples = true;

% Parse parameters
% for i=1:2:length(varargin)
%     switch(lower(varargin{i}))
%         case 'speedthresh'
%             speed_thresh = varargin{i+1};
%             if ~isa(speed_thresh, 'numeric') && speed_thresh <= 0
%                 error('Incorrect value for property ''SpeedThresh'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
%             end
%          case 'binsize'
%             binsize = varargin{i+1};
%             if ~isa(binsize, 'numeric')
%                 error('Incorrect value for property ''BinSize'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
%             end
%          case 'overlapfactor'
%             overlapFactor = varargin{i+1};
%             if ~isa(overlapFactor, 'numeric')
%                 error('Incorrect value for property ''overlapFactor'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
%             end
%          case 'savefigure'
%             savefig = varargin{i+1};
%             if savefig ~= 1 && savefig ~= 0
%                 error('Incorrect value for property ''SaveFigure'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
%             end
%          case 'plotexamples'
%             plotexamples = varargin{i+1};
%             if plotexamples ~= 1 && plotexamples ~= 0
%                 error('Incorrect value for property ''PlotExamples'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
%             end
%          case 'restrictsleep'
%             SleepTimeToRestrict = varargin{i+1};
%             if ~isa(SleepTimeToRestrict, 'numeric')
%                 error('Incorrect value for property ''RestrictSleep'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
%             end
%          case 'repeatoriginalpaper'
%             do_as_in_paper = varargin{i+1};
%             if do_as_in_paper ~= 1 && do_as_in_paper ~= 0
%                 error('Incorrect value for property ''RepeatOriginalPaper'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
%             end
%     end
% end
%% Load data
Dir = PathForExperimentsERC_Arsenii('UMazePAG');
Dir = RestrictPathForExperiment(Dir, 'nMice', nmice); %restrict on mice, which are used in this session (nmice)

%% Allocate memory

spikes = cell(1, length(Dir.path)); %ячейка (массив) размером 1 на длину соответствующей строки из load`a
behav = cell(1, length(Dir.path));
sleep = cell(1, length(Dir.path));
Epochs = cell(1, length(Dir.path));% Epoch - вектор, имеющий начало и конец временного интервала
UMazeEpoch = cell(1, length(Dir.path));

overlappairs = cell(1, length(Dir.path));
distantpairs = cell(1, length(Dir.path));

S_PC = cell(1, length(Dir.path));
QEpoched = cell(1, length(Dir.path));

rho_overlap = cell(1, length(states));
rho_distant = cell(1, length(states));

for iepoch = 1:length(states)
    rho_overlap{iepoch} = nan(1, 10e5);
    rho_distant{iepoch} = nan(1, 10e5);
end

for i = 1:length(Dir.path)
    spikes{i} = load([Dir.path{i}{1} 'SpikeData.mat'], 'S', 'PlaceCells', 'TT', 'RippleGroups'); %Загружаем данные спайков. load([путь, название файла], переменные); S - spiketrains, PlaceCells - indices of PC
    behav{i} = load([Dir.path{i}{1} 'behavResources.mat'], 'SessionEpoch', 'CleanVtsd', 'CleanAlignedXtsd', 'CleanAlignedYtsd', 'FreezeAccEpoch'); %Загружаем поведенческие данные. load([путь, название файла], переменные);
    if strcmp(Dir.name{i}, 'Mouse906') || strcmp(Dir.name{i}, 'Mouse977') || strcmp(Dir.name{i}, 'Mouse861') || strcmp(Dir.name{i}, 'Mouse979') %Mice with bad OB-based sleep scoring.
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_Accelero.mat'], 'SWSEpoch'); %SleepScoring is the process, which helps you distinguish the epoch.
    else
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_OBGamma.mat'], 'SWSEpoch'); %SleepScoring происходит либо по гамма-осцилляциям обонятельной луковицы, либо по ускорению 
    end
end    

%% Create epochs
for i = 1:length(Dir.path)
    [~, UMazeEpoch{i}, ConditioningEpoch, TaskEpoch, AfterConditioningEpoch] = ReturnMnemozyneEpochs(behav{i}.SessionEpoch, ...
        'Speed', behav{i}.CleanVtsd, 'SpeedThresh', speed_thresh);
    
    % Get freezing epoch during cond
    [~, ~, tempEpoch] = ReturnMnemozyneEpochs(behav{i}.SessionEpoch);
    ConditioningFreezingEpoch = and(behav{i}.FreezeAccEpoch, tempEpoch);
    
    % Sleep Epochs
    if isempty(SleepTimeToRestrict)
        PreSleepFinal = and(behav{i}.SessionEpoch.PreSleep, sleep{i}.SWSEpoch);
        PostSleepFinal = and(behav{i}.SessionEpoch.PostSleep, sleep{i}.SWSEpoch);
    else
        PreSleepFinal = RestrictToTime(and(behav{i}.SessionEpoch.PreSleep, sleep{i}.SWSEpoch), SleepTimeToRestrict);
        PostSleepFinal = RestrictToTime(and(behav{i}.SessionEpoch.PostSleep, sleep{i}.SWSEpoch), SleepTimeToRestrict);
    end
    
    % Pack everything in one cell
    Epochs{i} = {PreSleepFinal, TaskEpoch, ConditioningEpoch, ConditioningFreezingEpoch, PostSleepFinal, AfterConditioningEpoch};
end

%% Create binned FR vectors. ZScore - нормализация. В этом случае по времени. Она нужна, поскольку firing rate нейронов разный. Получаем отклонение каждого нейрона от среднего firing rate по этому нейрону. Из каждого значения вычитаем среднее и делим на стандартное отклонение. Для каждой эпохи создаем Spike histograms
for i=1:length(Dir.path)
    QEpoched{i} = cell(1, length(states));
%     id_ripples{i} = PyramLayerSorting(spikes{i}.S, spikes{i}.RippleGroups, spikes{i}.TT);
%     pyr_cells = spikes{i}.S(id_ripples{i});
        
    if do_as_in_paper
        Q = MakeQfromS(spikes{i}.S, 10); % 1 ms
        
        for iepoch = 1:length(states)
            QEpoched{i}{iepoch} = full(Data(Restrict(Q, Epochs{i}{iepoch})));
        end       
    else
        Q = MakeQfromS(spikes{i}.S, binsize);
        
        for iepoch = 1:length(states)
            QEpoched{i}{iepoch} = zscore(full(Data(Restrict(Q, Epochs{i}{iepoch}))));
        end   
    end
    clear Q

end

%% Calculate ovelap for each place field
for i = 1:length(Dir.path)
    %Get the stats for each cell
    if isfield(spikes{i}.PlaceCells, 'idx') %Есть ли у мыши клетки места. То есть есть ли у нее в PlaceFields 'idx'
        if length(spikes{i}.PlaceCells.idx) > PCnum_thresh %Take only those mice which have more than PCnum_thresh PC
            S_PC{i} = spikes{i}.S(spikes{i}.PlaceCells.idx); %Из spikes.S берем только индексы клеток места
            [overlappairs{i}, distantpairs{i}] = FindOverlappingPlaceFields(S_PC{i}, behav{i}.CleanAlignedXtsd, behav{i}.CleanAlignedYtsd, UMazeEpoch{i}, overlapFactor);
        else
            overlappairs{i}={};
            distantpairs{i}={};
            fprintf(['Error. length(spikes{i}.PlaceCells.idx) <= PCnum_thresh in mouse ' num2str(nmice(i)) '\n']);
        end
    else
         overlappairs{i}={};
         distantpairs{i}={};
         fprintf('Error #1 if')
    end
end 

%% Look at each pair overlapping
% if plotexamples
%     for imouse = 1:length(Dir.path)
%         PlotCorrelationExamples(Dir.name{imouse}, spikes{imouse}, S_PC{imouse}, overlappairs{imouse},...
%             Epochs{imouse}{1}, Epochs{imouse}{2}, Epochs{imouse}{5});
%     end
% end

%% Calculate cross-correlation between spike trains

cnt_o = 1;
cnt_d = 1;

for i = 1:length(Dir.path)
    
    if do_as_in_paper
        % Perform calculations - overlapping cells
        if ~isempty(overlappairs{i})
            for j=1:length(overlappairs{i})
                pair = overlappairs{i}{j};
                for iepoch = 1:length(states)
                    rho_overlap{iepoch}(cnt_o) = mean(crosscorr(QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                        QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(2))), 'NumLags', 100));
%                     rho_overlap{iepoch}(cnt_o) = mean(CrossCorrDB(Data(Restrict(S_PC{i}{pair(1)}, Epochs{i}{iepoch})),...
%                         Data(Restrict(S_PC{i}{pair(2)}, Epochs{i}{iepoch})), 2, 100));
                end
                cnt_o = cnt_o+1;
            end
        end
        clear pair
        % Perform calculations - non-overlapping cells
        if ~isempty(distantpairs{i})
            for j=1:length(distantpairs{i})
                pair = distantpairs{i}{j};
                for iepoch = 1:length(states)
                    rho_distant{iepoch}(cnt_d) = mean(crosscorr(QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                        QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(2))), 'NumLags', 100));
                end
                cnt_d=cnt_d+1;
            end
        end
        
    else % do_as_in_paper
        if ~isempty(overlappairs{i}) %Calculate ovelapping cells
            for j = 1:length(overlappairs{i})
                pair = overlappairs{i}{j};
                for iepoch = 1:length(states)
                    temp = corrcoef(QEpoched{i}{iepoch}(:, spikes{i}.PlaceCells.idx(pair(1))), QEpoched{i}{iepoch}(:, spikes{i}.PlaceCells.idx(pair(2))));
                    rho_overlap{iepoch}(cnt_o) = temp(1,2);
                    clear temp;
                end
                cnt_o = cnt_o + 1;
            end        
        end
        clear pair;
    
        if ~isempty(distantpairs{i}) %Calculate non-ovelapping cells
            for j = 1:length(distantpairs{i})
                pair = distantpairs{i}{j};
                for iepoch = 1:length(states)
                    temp = corrcoef(QEpoched{i}{iepoch}(:, spikes{i}.PlaceCells.idx(pair(1))), QEpoched{i}{iepoch}(:, spikes{i}.PlaceCells.idx(pair(2))));
                    rho_distant{iepoch}(cnt_d) = temp(1,2);
                    clear temp;
                end
                cnt_d = cnt_d + 1;
            end        
        end
    end
end

% Remove missing values (PostSleep is a reference)
idx_toremove_overlap = find(isnan(rho_overlap{5}),1);
idx_toremove_distant = find(isnan(rho_distant{5}),1);
for iepoch = 1:length(states)
    rho_overlap{iepoch}(idx_toremove_overlap:end) = [];
    rho_distant{iepoch}(idx_toremove_distant:end) = [];
end

%% Figures
% Short figure 1
foldertosave = ChooseFolderForFigures_AG('ReactReplay');
f1 = figure('units', 'normalized', 'outerposition', [0 0 0.4 0.7]);
cols_ssd = cols;
while length(cols_ssd) < 6
    cols_ssd{length(cols_ssd)+1} = cols{1};
    cols_ssd{length(cols_ssd)+1} = cols{2};
end

short_stacked_data = {rho_overlap{1}, rho_distant{1}, rho_overlap{3}, rho_distant{3}, rho_overlap{5}, rho_distant{5}};
ssp = MakeBoxPlot_DB(short_stacked_data, cols_ssd, [1 2 4 5 7 8], [], 0);
set(gca,'Xtick',1.5:3:7.5,'XtickLabel',{'PreSleep', 'Conditioning', 'PostSleep'}, 'FontSize', 15)
ylabel('Corellation', 'FontSize', 15)
title('Pairwise correlations', 'FontSize', 15);
legend([ssp{1}.handles.box, ssp{2}.handles.box], {'Overlapping PCs', 'Non-overlapping PCs'});

pairs_toCompare = {[1 3], [1 5], [3 5]};
% set(pairs_toCompare{1}, 'Xtick', [1 3]);


% p = DoWilcoxonOnArray(short_stacked_data, pairs_toCompare);
% for ip = 1:length(p)
%     if p(ip) <= 0.05
%         sigstar_DB(pairs_toCompare(ip),p(ip),0,'LineWigth',16,'StarSize',24);
%     end
% end
% 
% pairs_toCompare = {[2 4], [2 6], [4 6]};
% p = DoWilcoxonOnArray(short_stacked_data, pairs_toCompare);
% for ip = 1:length(p)
%     if p(ip) <= 0.05
%         sigstar_DB(pairs_toCompare(ip),p(ip),0,'LineWigth',16,'StarSize',24);
%     end
% end
% 
% makepretty
% % Save the figure
% if savefig
%     saveas(f1,[foldertosave '/Pairwise_Small.fig']);
%     saveFigure(f1,'Pairwise_Small', foldertosave);
% end

%% Final figures for Karim
% Reorganize data a bit (PreSleep-PostSleep-Maze-CondMov-CondFreeze-PostTests)
largefig_data{1} = {rho_overlap{1}, rho_overlap{5}, rho_overlap{2}, rho_overlap{3}, rho_overlap{4}, rho_overlap{6}};
largefig_data{2} = {rho_distant{1}, rho_distant{5}, rho_distant{2}, rho_distant{3}, rho_distant{4}, rho_distant{6}};
pairs_toCompare = {[1 2], [1 3], [2 3]};

%% Big figure
f2 = figure('units', 'normalized', 'outerposition', [0.1 0.1 0.83 0.7]);
ax{1} = axes('position', [0.05 0.05 0.43 0.9]);
ax{2} = axes('position', [0.54 0.05 0.43 0.9]);
titles = {'Overlapping place cells', 'Non-Overlapping place cells'};
for icolor = 1:length(cols)
    cols_lf{icolor} = repelem(cols(icolor), length(states));
end
% Plot
for iaxis = 1:length(ax)
    axes(ax{iaxis});
    MakeBoxPlot_DB(largefig_data{iaxis}, cols_lf{iaxis}, 1:length(states), [], 0);
    set(gca,'Xtick', 1:6, 'XtickLabel',{'PreSleep', 'PostSleep', 'NoStimTask', 'StimRun', 'StimFreeze', 'PostTest'})
    ylabel('Correlation', 'FontSize', 18)
    title(titles{iaxis});
    p = DoWilcoxonOnArray(largefig_data{iaxis}, pairs_toCompare);
    for ip = 1:length(p)
        if p(ip) <= 0.05
            sigstar_DB(pairs_toCompare(ip),p(ip),0,'LineWigth',16,'StarSize',24);
        end
    end
%     makepretty
end
ylims_out = SelectYlim(f2);
for iaxis = 1:length(ax)
    axes(ax{iaxis});
    ylim(ylims_out);
end
if savefig
    saveas(f2,[foldertosave '/Pairwise_final_big.fig']);
    saveFigure(f2,'Pairwise_final_big', foldertosave);
end

%% Short figure 2
shortfig_data{1} = {rho_overlap{1}, rho_overlap{2}, rho_overlap{5}};
shortfig_data{2} = {rho_distant{1}, rho_distant{2}, rho_distant{5}};

f3 = figure('units', 'normalized', 'outerposition', [0.1 0.15 0.5 0.7]);
ax{1} = axes('position', [0.05 0.05 0.43 0.9]);
ax{2} = axes('position', [0.54 0.05 0.43 0.9]);
for icolor = 1:length(cols)
    cols_sf{icolor} = repelem(cols(icolor), 3);
end
for iaxis = 1:length(ax)
    axes(ax{iaxis});
    MakeBoxPlot_DB(shortfig_data{iaxis}, cols_sf{iaxis}, 1:3, [], 0);
    set(gca,'Xtick', 1:3, 'XtickLabel',{'PreSleep', 'Full Task', 'PostSleep'})
    ylabel('Correlation', 'FontSize', 16)
    title(titles{iaxis});
    p = DoWilcoxonOnArray(shortfig_data{iaxis}, pairs_toCompare);
    for ip = 1:length(p)
        if p(ip) <= 0.05
            sigstar_DB(pairs_toCompare(ip),p(ip),0,'LineWigth',16,'StarSize',24);
        end
    end
    makepretty
end
ylims_out = SelectYlim(f3);
for iaxis = 1:length(ax)
    axes(ax{iaxis});
    ylim(ylims_out);
end
if savefig
    saveas(f3,[foldertosave '/Pairwise_final_small.fig']);
    saveFigure(f3,'Pairwise_final_small', foldertosave);
end

%% Only overlapping cells
f4 = figure('units', 'normalized', 'outerposition', [0.3 0.2 0.285 0.7]);
MakeBoxPlot_DB(shortfig_data{1}, cols_sf{1}, 1:3, [], 0);
set(gca,'Xtick', 1:3, 'XtickLabel',{'PreSleep', 'Full Task', 'PostSleep'})
ylabel('Correlation', 'FontSize', 18)
p = DoWilcoxonOnArray(shortfig_data{1}, pairs_toCompare);
for ip = 1:length(p)
    if p(ip) <= 0.05
        sigstar_DB(pairs_toCompare(ip),p(ip),0,'LineWigth',16,'StarSize',24);
    end
end
% makepretty
if savefig
    saveas(f4,[foldertosave '/Pairwise_fullfinal_pooled.fig']);
    saveFigure(f4,'Pairwise_fullfinal_pooled', foldertosave);
end


%% Auxiliary functions

function PlotCorrelationExamples(MouseName, spikes, S_PC, overlappairs, PreSleepEpoch, TaskEpoch, PostSleepEpoch, varargin)
    
% Optional parameters
binsize = 4;
windowsize = 400;
% Parse parameters
for k=1:2:length(varargin)
    switch(lower(varargin{k}))
         case 'binsize'
            binsize = varargin{k+1};
            if ~isa(binsize, 'numeric')
                error('Incorrect value for property ''BinSize'' (type ''help PlotCorrelationExamples'' for details).');
            end
         case 'windowsize'
            windowsize = varargin{k+1};
            if ~isa(windowsize, 'numeric')
                error('Incorrect value for property ''WindowSize'' (type ''help PlotCorrelationExamples'' for details).');
            end
    end
end

nbins = windowsize*2/binsize;

% Do the job and plot the figure

if ~isempty(overlappairs)
    CPre = cell(length(overlappairs),1);
    CTask = cell(length(overlappairs),1);
    CPost = cell(length(overlappairs),1);
    for ipair=1:length(overlappairs)
        [CPre{ipair},B] = CrossCorrDB(Data(Restrict(S_PC{overlappairs{ipair}(1)},PreSleepEpoch)),...
            Data(Restrict(S_PC{overlappairs{ipair}(2)},PreSleepEpoch)),binsize,nbins);
        CTask{ipair} = CrossCorrDB(Data(Restrict(S_PC{overlappairs{ipair}(1)},TaskEpoch)),...
            Data(Restrict(S_PC{overlappairs{ipair}(2)},TaskEpoch)),binsize,nbins);
        CPost{ipair} = CrossCorrDB(Data(Restrict(S_PC{overlappairs{ipair}(1)},PostSleepEpoch)),...
            Data(Restrict(S_PC{overlappairs{ipair}(2)},PostSleepEpoch)),binsize,nbins);
    end
    
    % Figure
    if length(overlappairs) < 15
        fff = figure('units', 'normalized', 'outerposition', [0 0.55 0.7 0.4], 'Name', MouseName);
        title(MouseName);
        tabs = arrayfun(@(x) uitab('Title', ['Cl' num2str(spikes.PlaceCells.idx(overlappairs{x}(1)))...
            ' vs Cl' num2str(spikes.PlaceCells.idx(overlappairs{x}(2)))]), 1:length(overlappairs));
        ax = arrayfun(@(tab) axes(tab), tabs);
        arrayfun(@(k) title(ax(k), ['Cl' num2str(spikes.PlaceCells.idx(overlappairs{k}(1)))...
            ' vs Cl' num2str(spikes.PlaceCells.idx(overlappairs{k}(2)))], 'FontSize', 16), 1:length(overlappairs));
        arrayfun(@(k) set(ax(k), 'NextPlot','add', 'FontSize', 16, 'FontWeight', 'bold'), 1:length(overlappairs));
        arrayfun(@(k) plot(ax(k), B, CPre{k}, 'k', 'LineWidth', 2),...
            1:length(overlappairs));
        arrayfun(@(k) plot(ax(k), B, CTask{k}, 'm', 'LineWidth', 2),...
            1:length(overlappairs));
        arrayfun(@(k) plot(ax(k), B, CPost{k}, 'b', 'LineWidth', 2),...
            1:length(overlappairs));
        arrayfun(@(k) line(ax(k),[0 0], ylim, 'Color', 'r', 'LineWidth',3), 1:length(overlappairs));
        arrayfun(@(k) ylabel(ax(k), 'Count'), 1:length(overlappairs));
        arrayfun(@(k) xlabel(ax(k), 'Time (ms)'), 1:length(overlappairs));
        arrayfun(@(k) legend(ax(k), {'Pre', 'Task', 'Post'}), 1:length(overlappairs));
    else
        for numfig = 1:ceil(length(overlappairs)/15)
            try
                rr = ((numfig-1)*15)+1:((numfig-1)*15)+15;
                temp_o = overlappairs(rr);
            catch
                rr = ((numfig-1)*15)+1:length(overlappairs);
                temp_o = overlappairs(rr);
            end
            fff1(numfig) = figure('units', 'normalized', 'outerposition', [0 0.55 0.7 0.4], 'Name', MouseName);
            title(MouseName);
            tabs = arrayfun(@(x) uitab('Title', ['Cl' num2str(spikes.PlaceCells.idx(temp_o{x}(1)))...
                ' vs Cl' num2str(spikes.PlaceCells.idx(temp_o{x}(2)))]), 1:length(temp_o));
            ax = arrayfun(@(tab) axes(tab), tabs);
            arrayfun(@(k) title(ax(k), ['Cl' num2str(spikes.PlaceCells.idx(temp_o{k}(1)))...
                ' vs Cl' num2str(spikes.PlaceCells.idx(temp_o{k}(2)))], 'FontSize', 16), 1:length(temp_o));
            arrayfun(@(k) set(ax(k), 'NextPlot','add', 'FontSize', 16, 'FontWeight', 'bold'), 1:length(temp_o));
            arrayfun(@(k) plot(ax(k+1-rr(1)), B, CPre{k}, 'k', 'LineWidth', 2),...
                rr);
            arrayfun(@(k) plot(ax(k+1-rr(1)), B, CTask{k}, 'm', 'LineWidth', 2),...
                rr);
            arrayfun(@(k) plot(ax(k+1-rr(1)), B, CPost{k}, 'b', 'LineWidth', 2),...
                rr);
            arrayfun(@(k) line(ax(k),[0 0], ylim, 'Color', 'r', 'LineWidth',3), 1:length(temp_o));
            arrayfun(@(k) ylabel(ax(k), 'Count'), 1:length(temp_o));
            arrayfun(@(k) xlabel(ax(k), 'Time (ms)'), 1:length(temp_o));
            arrayfun(@(k) legend(ax(k), {'Pre', 'Task', 'Post'}), 1:length(temp_o));
        end
    end
end
end