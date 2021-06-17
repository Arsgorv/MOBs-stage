% This function plots spike trains imagesc

Q_S = cell(1, length(Dir.path));
Data_S = cell(1, length(Dir.path));

cnt = 0;
for ipic = 1:length(Dir.path)
    if cnt == 0
        figure('units', 'normalized', 'outerposition', [0 0 1 0.5]);
    end
        cnt = cnt + 1;
        subplot(1, 3, cnt)
               
        % Load data
        BR_S = load([Dir.path{ipic}{1} 'behavResources.mat'], 'SessionEpoch');
        
        % Perform spike binning and plotting
        Q_S{ipic} = MakeQfromS(s{ipic}.S, 0.1e4);
%         timeQ = Range(Q_S{ipic}, 's');
        Data_S{ipic} = zscore(full(Data(Q_S{ipic})));
        
%         f1_subplot = imagesc(timeQ, 1:size(Data_S{ipic}'), Data_S{ipic}');
        f1_subplot = imagesc(Data_S{ipic}');
        caxis([-1 3]);
        ylabel('number of neurons');
        xlabel('time (number of bins)');
        title([num2str(Dir.name{ipic})], 'FontSize', 14);
        
        % Separate sessions on the plot.
        % White lines - PreSleep Epoch,
        % Red lines - PostSleep Epoch
                
        PreSleep_start_time_point = Start(BR_S.SessionEpoch.PreSleep)/0.1e4;
        PreSleep_end_time_point = End(BR_S.SessionEpoch.PreSleep)/0.1e4;
        PostSleep_start_time_point = Start(BR_S.SessionEpoch.PostSleep)/0.1e4;
        PostSleep_end_time_point = End(BR_S.SessionEpoch.PostSleep)/0.1e4;
        hold on
        line([PreSleep_start_time_point PreSleep_start_time_point], ylim, 'Color', [1 1 1], 'LineWidth', 2);
        line([PreSleep_end_time_point PreSleep_end_time_point], ylim, 'Color', [1 1 1], 'LineWidth', 2);
        line([PostSleep_start_time_point PostSleep_start_time_point], ylim, 'Color', [1 0 0], 'LineWidth', 2);
        line([PostSleep_end_time_point PostSleep_end_time_point], ylim, 'Color', [1 0 0], 'LineWidth', 2);
        
    if cnt == 3
        cnt = 0;
    end  
end