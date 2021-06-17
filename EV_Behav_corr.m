function [f1,f2] = EV_Behav_corr(nMice, varargin)
%
% This function calculates Spearman correlation between explained variance and behavioral parameters, such as:
%
%     * PostSleep - PreSleep difference of the first entry to shock or safe zone time.
%     * PostSleep - PreSleep difference of the number of entries to shock or safe zones.
%     * PostSleep - PreSleep occupation difference of shock or safe zones.
%     * PostSleep - PreSleep speed difference of shock or safe zones.
%
% And plots it as scatter plots with least-squares line.
%
% Take care of the organization of input EV and behavior files. They should respect the order of the mice (nMice)
% As for EV, it should be structured as EV{NREM/REM, 1}{Epoch, 1}(nMouse).
% To calculate EV please use ExplainedVariance_new_AG.m
%
% As for behavior:
%
% val(:, :, 1) = 
% mouse 1    shock_pre - safe_pre    shock_post - safe_post
% mouse 2    shock_pre - safe_pre    shock_post - safe_post
% ...
%
% However, originally, I correlate EV/REV only with POST-PRE differences:
% val(:, :, 2) =
% mouse 1    post_shock - pre_shock    post_safe - pre_safe
% mouse 2    post_shock - pre_shock    post_safe - pre_safe
%
% Coded by Arsenii Goriachenkov, MOBS team, Paris,
% 29/03/2021
%
% github.com/Arsgorv

%% Temporary parameters
% AddMyPaths_Arsenii;

% nMice = [905 906 911 994 1161 1162 1168 1182 1186]; 
nMice = [905 906 911 994 1161 1162 1168];

%% Load data
% foldertoload = ChooseFolderForFigures_DB('Data');

% If you want to make analysis with cells of pyramidal layer and interneurons, load this:
% EV_Res = load(fullfile(foldertoload, 'EV_res.mat'));

% If you want to make analysis with cells of just pyramidal layer, load this:
% EV_Res = load('EV_res_pyr.mat');
EV_Res = load('EV_res_pyr_inter.mat');


% Load behavior parameters
% Behav_ratios = load(fullfile(foldertoload, 'Behav_ratios.mat'));
Behav_ratios = load('Behav_ratios.mat');

% To check differences between EV and REV
% EV_Res = load('EV_REV_diff.mat');
% EV_Res.EV_REV_diff{1} = EV_Res.EV_REV_diff{1}';
% EV_Res.EV_REV_diff{2} = EV_Res.EV_REV_diff{2}';
% EV_Res.REV_EV_diff{1} = EV_Res.REV_EV_diff{1}';
% EV_Res.REV_EV_diff{2} = EV_Res.REV_EV_diff{2}';

% To check ratio between EV and REV
% EV_Res = load('EV_REV_ratio.mat');
% EV_Res.EV_REV_ratio{1} = EV_Res.EV_REV_ratio{1}';
% EV_Res.EV_REV_ratio{2} = EV_Res.EV_REV_ratio{2}';
% EV_Res.REV_EV_ratio{1} = EV_Res.REV_EV_ratio{1}';
% EV_Res.REV_EV_ratio{2} = EV_Res.REV_EV_ratio{2}';

%% Choose epochs
list_epochs_wake = {'Explo', 'CondMov', 'CondFreeze', 'FullTask', 'RipplesEpoch', 'TestPost'};
[idx_wake_epoch,~] = listdlg('PromptString', {'Choose wakefullness epoch to build correlation matrices'},'ListString',list_epochs_wake, 'SelectionMode', 'single');
wake_epoch_name = list_epochs_wake{idx_wake_epoch};

list_epochs_sleep = {'NREM', 'REM'};
[idx_sleep_epoch,~] = listdlg('PromptString', {'Choose sleep type to build correlation matrices'},'ListString',list_epochs_sleep, 'SelectionMode', 'single');
sleep_epoch_name = list_epochs_sleep{idx_sleep_epoch};

%% Additional parameters
Zones = {'shock', 'safe'};
Behavior = {Behav_ratios.diff_fentry, Behav_ratios.diff_nentries, Behav_ratios.diff_occup, Behav_ratios.diff_speed};
Correlations = nan(length(Behavior), 2, length(Zones), 2, length(list_epochs_sleep)); %([fentry, nentries, occup, speed], [rho/p-value], [shock/safe zone], [EV/REV], [NREM/REM])
Names_for_figs = {' first entries time to ', ' number of entries to ', ' occupancy of ', ' average speed in ' };

Variance = {EV_Res.EV, EV_Res.REV};
Variance_names = {'EV', 'REV'};

% if you want to use differences or ratios of EV and REV instead of their original values, use this:
% Variance = {EV_Res.EV_REV_diff, EV_Res.REV_EV_diff,};
% Variance = {EV_Res.EV_REV_ratio, EV_Res.REV_EV_ratio};

% Variance_names = {'EV/REV difference', 'REV/EV difference'};
% Variance_names = {'EV/REV ratio', 'REV/EV ratio'};

colors = linspace(1,length(nMice),length(nMice));

%% Calculate correlations
for ibehav = 1:length(Behavior) % 1 - fentry, 2 - nentires, 3 - occup, 4 - speed
    for izone = 1:2 % 1 - shock zone, 2 - safe zone
        for ivariance = 1:2 % 1 - EV, 2 - REV
            [rho, p] = corr(Behavior{ibehav}(:, izone, 2), (Variance{ivariance}{idx_sleep_epoch}{idx_wake_epoch})', 'Type', 'Spearman', 'rows', 'complete');
            
            Correlations(ibehav, 1, izone, ivariance, idx_sleep_epoch) = rho;
            Correlations(ibehav, 2, izone, ivariance, idx_sleep_epoch) = p;
        end
    end
end

%% Plot figures
for ibehav = 1:2
    for ivariance = 1:length(Variance)
        f{ivariance} = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);
        top_left_axes = axes('position', [0.08 0.58 0.35 0.35]);
        top_right_axes = axes('position', [0.58 0.58 0.35 0.35]);
        bottom_left_axes = axes('position', [0.08 0.1 0.35 0.35]);
        bottom_right_axes = axes('position', [0.58 0.1 0.35 0.35]);
        
        Axes = [top_left_axes, top_right_axes, bottom_left_axes, bottom_right_axes];
        
        for ientry = 1:2
            for izone = 1:2
                if ientry == 1
                    axes(Axes(izone));
                elseif ientry == 2
                    axes(Axes(izone+2));
                end
                if ibehav == 1
                    temp = ientry;
                elseif ibehav == 2
                    temp = ientry + 2;
                elseif ibehav == 3
                    temp = ientry + 3;
                end
                scatter(Variance{ivariance}{idx_sleep_epoch}{idx_wake_epoch}, Behavior{temp}(:, izone, 2), 75, 'k', 'filled')
                title(['Post - Pre differences of the' Names_for_figs{temp} ' the ' Zones{izone} ' zone']);
                
                t1(izone) = text(.4, .9, ['rho = ' num2str(round(Correlations(temp, 1, izone, ivariance, idx_sleep_epoch), 3))], 'sc', 'FontSize', 13);
                t2(izone) = text(.4, .8, ['p = ' num2str(round(Correlations(temp, 2, izone, ivariance, idx_sleep_epoch), 3))], 'sc', 'FontSize', 13);
                l1 = lsline;
                l1.Color = 'k';
                l1.LineWidth = 1;
                xlabel([Variance_names{ivariance} ' ' wake_epoch_name '. ' sleep_epoch_name]);
                ylabel('Post - Pre differences');
                set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
                set(gca, 'LineWidth', 1);
            end
        end
    end
end

end