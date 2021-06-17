% This function plots matrices of correlations for Pre-, Task and Post for CorrM
% You can choose the type of sleep: either NREM or REM and the wakefullness epoch
%
% By Arsenii Goriachenkov, MOBS team, Paris,
% 04/05/2021
% github.com/arsgorv

list_epochs_wake = {'Explo', 'CondMov', 'CondFreeze', 'FullTask', 'TestPost'};
[idx_wake_epoch,~] = listdlg('PromptString', {'Choose wakefullness epoch to build correlation matrices'},'ListString',list_epochs_wake, 'SelectionMode', 'single');
wake_epoch_name = list_epochs_wake{idx_wake_epoch};

list_epochs_sleep = {'NREM', 'REM'};
[idx_sleep_epoch,~] = listdlg('PromptString', {'Choose sleep type to build correlation matrices'},'ListString',list_epochs_sleep, 'SelectionMode', 'single');
sleep_epoch_name = list_epochs_sleep{idx_sleep_epoch};

cnt = 0;
for i = 1:length(Dir.path)
    if cnt == 0
        figure('units', 'normalized', 'outerposition', [0 0 0.9 0.5]);
    end
    cnt = cnt + 1;
  
    % Plot NREM correlations. Full Task
    subplot(1, 3, cnt)
    imagesc(CorrM{idx_sleep_epoch, 1}{1, idx_wake_epoch}(i).pre)
    caxis([-0, .1]);
    ylabel('neuron ¹');
    xlabel('neuron ¹');
    title([num2str(Dir.name{i}) ' PreSleep. ' sleep_epoch_name ' correlation'], 'FontSize', 14);
    
    subplot(1, 3, cnt + 1)
    imagesc(CorrM{idx_sleep_epoch, 1}{1, idx_wake_epoch}(i).task)
    caxis([-0, .1]);
    ylabel('neuron ¹');
    xlabel('neuron ¹');
    title([num2str(Dir.name{i}),' ' wake_epoch_name '. ' sleep_epoch_name ' correlation'], 'FontSize', 14);
    
    subplot(1, 3, cnt + 2)
    imagesc(CorrM{idx_sleep_epoch, 1}{1, idx_wake_epoch}(i).post)
    caxis([-0, .1]);
    ylabel('neuron ¹');
    xlabel('neuron ¹');
    title([num2str(Dir.name{i}) ' PostSleep. ' sleep_epoch_name ' correlation'], 'FontSize', 14);

    cnt = 0;
end