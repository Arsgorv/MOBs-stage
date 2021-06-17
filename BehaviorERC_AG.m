function fh = BehaviorERC_AG(nMice, varargin)
%
% This function plots 4 behavioral metrics of the Shock Zone in the UMaze
% for Pre- and Post-Tests:
%  - Occupancy of the zone
%  - Number of entries into the zone
%  - Time to enter the zone for the first time
%  - Average speed in the zone
%
% INPUT
%
%     nMice      numbers of mice to include in the analysis
%     IsSave     (optional) true or 1 if to save the figure (default=false)
%     SownPoints (optional) 1 if you want to show individual points (default=0)
% 
%  OUTPUT
%
%     nothing
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 05/11/2020 based on the script I wrote in 2018
% github.com/bryzgalovdm
% 
% UPDATED:
% Added:
%           - [shock_pre/safe_pre, shock_post/safe_post]; [shock_pre/shock_post, safe_pre/safe_post]; There are both ratios and differences for each parameter
%           - Safe-zone occupancy, Safe-zone speed. 
%           - Plots for every ratio or every difference.
%           - FT = [0/1] - you can choose to look at only 2 first trials (ex. TestPost1, TestPost2) or to all four trials
%           - Calculate speed distribution
%           - Plot mice trajectories in Pre/Post tests
%
% By Arsenii Goriachenkov, MOBS team, Paris,
% 23/03/2021
% github.com/arsgorv

%% Temporary parameters
AddMyPaths_Arsenii;

% nMice = [797 798 828 861 882 905 906 911 912 977 994 1117 1124 1161 1162 1168 1182 1186];
% nMice = [905 906 911 994 1161 1162 1168 1182 1186]; 
% nMice = [905 906 911 994 1161 1162 1168]; 
% nMice = [1182 1186]; 

%% Parameters
FigName = 'UMazePAG_Behavior';
IsSave = false;
SP = 1; %show points
% Do you want to save the result?
save_res = 1; % Put '1' if yes.
% foldertosave = ChooseFolderForFigures_DB('Data');
foldertosave = 'E:\ERC_data';

% foldertosave = ChooseFolderForFigures_DB('Data');

FT = 0; % 1 - first two trials (ex. TestPost1, TestPost2). 0 - all trials (TestPost1, TestPost2, TestPost3, TestPost4)
if FT == 0
    num_sess = 4;
elseif FT == 1
    num_sess = 2;
end

%% Optional Parameters
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            IsSave = varargin{i+1};
            if IsSave ~= 1 && IsSave ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help PreTestCharacteristics'' for details).');
            end
        case 'showpoints'
            SP = varargin{i+1};
            if SP ~= 1 && SP ~= 0
                error('Incorrect value for property ''ShowPoints'' (type ''help PreTestCharacteristics'' for details).');
            end
    end
end

%% Get data
% Dir = PathForExperimentsERC('UMazePAG');
Dir = PathForExperimentsERC_Arsenii('UMazePAG');

Dir = RestrictPathForExperiment(Dir,'nMice', nMice);

a = cell(length(Dir.path),1); % Each "a" cell has behavResources of a certain mouse
for i = 1:length(Dir.path)
    a{i} = load([Dir.path{i}{1} '/behavResources.mat'], 'behavResources');
end
clear i

%% Find indices of PreTests and PostTest session in the structure
id_Pre = cell(1,length(a));
id_Post = cell(1,length(a));

for i=1:length(a) %find index of a certain epoch in behavResources
    id_Pre{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPre');
    id_Post{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPost');
    
    if FT == 0
        Test_iterr = length(id_Pre{i});
    elseif FT == 1
        Test_iterr = 2;
    end
end
clear i

%% Calculate average occupancy in shock and save zones
% Calculate occupancy of shock zone � safe zone for each mouse for each of TestPre/TestPost 

occup_shock = nan(length(Dir.path), num_sess, 2); % 4 tests (TestPre1, TestPre2, TestPre3, TestPre4 for example), Pre and Post
occup_safe = nan(length(Dir.path), num_sess, 2); % 4 tests, Pre and Post

for i=1:length(a) %Mouse itterator
    for k=1:Test_iterr %Test itterator (4 trials)
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Pre{i}(k))); %Output - occup - vector of seven cells, each of which - occupancy time of each epoch.
        occup_shock(i,k,1) = temp(1);
        occup_safe(i,k,1) = temp(2);
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Post{i}(k)));
        occup_shock(i,k,2) = temp(1);
        occup_safe(i,k,2) = temp(2);
    end
end

% Mean (of tests, not of mice)
occup_mean = nan(length(Dir.path), 2, 2); 
occup_std = nan(length(Dir.path), 2, 2);

for izone = 1:2 % 1st val codes preTest, 2nd val codes postTest; izone = 1 - Pre test, izone = 2 - Post test
    occup_mean(:, 1, izone) = mean(squeeze(occup_shock(:,:,izone)), 2); %Filling shock
    occup_mean(:, 2, izone) = mean(squeeze(occup_safe(:,:,izone)), 2); %Filling safe
    
    occup_std(:, 1, izone) = std(squeeze(occup_shock(:,:,izone)),0,2); %Filling shock
    occup_std(:, 2, izone) = std(squeeze(occup_safe(:,:,izone)),0,2); %Filling safe
end

%% calculate ratios of occupancy
% structure of ratio_occup:
% val(:, :, 1) = 
%             Pre       Post
% mouse 1 shock/safe shock/safe
% mouse 2 shock/safe shock/safe
% ...
% 
% val(:, :, 2) = 
%                 Shock               Safe
% mouse 1 shock_post/shock_pre safe_post/safe_pre
% mouse 2 shock_post/shock_pre safe_post/safe_pre
% ...

ratio_occup = nan(length(Dir.path), 2, 2);
diff_occup = nan(length(Dir.path), 2, 2);
for i=1:length(Dir.path)  % fix mouse (row)
    for j = 1:2 % fix either Pre/Post val, or Shock/Safe val
        for k = 1:2 % fix column of Pre/Post in 1st val or Shock/Safe in 2nd val
            if j == 1 %Val 1: columns - Pre, Post. Cell - shock/safe
                ratio_occup(i, k, j) = occup_mean(i, 1, k)/occup_mean(i, 2, k); 
                diff_occup(i, k, j) = occup_mean(i, 1, k)-occup_mean(i, 2, k);
            elseif j == 2 %Val 2: columns - Shock, Safe. Cell - Post/Pre
                ratio_occup(i, k, j) = occup_mean(i, k, 2)/occup_mean(i, k, 1);
                diff_occup(i, k, j) = occup_mean(i, k, 2)-occup_mean(i, k, 1);
            end
        end
    end
end

%% Calculate the 'first entry to the shock/safe zone' time

EntryTime_shock = nan(length(Dir.path), num_sess, 2); % 4 tests, Pre and Post
EntryTime_safe = nan(length(Dir.path), num_sess, 2); % 4 tests, Pre and Post

for i = 1:length(a)
    for k=1:Test_iterr % Iterator for tests
        temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Pre{i}(k)), 240);
        for j = 1:length(temp) % Substitute 240 for NaN. 
            if temp(j) == 240
                temp(j) = NaN;
            end
        end
        EntryTime_shock(i,k,1) = temp(1);
        EntryTime_safe(i,k,1) = temp(2);
        temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Post{i}(k)), 240);
        for j = 1:length(temp)
            if temp(j) == 240 % Substitute 240 for NaN. 
                temp(j) = NaN;
            end
        end
        EntryTime_shock(i,k,2) = temp(1);
        EntryTime_safe(i,k,2) = temp(2);
    end
end

EntryTime_mean = nan(length(Dir.path), 2, 2); 
EntryTime_std = nan(length(Dir.path), 2, 2); 

for izone = 1:2 % 1st val codes preTest, 2nd val codes postTest
    EntryTime_mean(:, 1, izone) = nanmean(squeeze(EntryTime_shock(:,:,izone)),2);
    EntryTime_mean(:, 2, izone) = nanmean(squeeze(EntryTime_safe(:,:,izone)),2);
    
    EntryTime_std(:, 1, izone) = std(squeeze(EntryTime_shock(:,:,izone)),0,2);
    EntryTime_std(:, 2, izone) = std(squeeze(EntryTime_safe(:,:,izone)),0,2);
end

%% Calculate ratios of 'first entry to the shock/safe zone'
% structure of ratio_fentry:
% val(:, :, 1) = 
%             Pre       Post
% mouse 1 shock/safe shock/safe
% mouse 2 shock/safe shock/safe
% ...
% 
% val(:, :, 2) = 
%                 Shock               Safe
% mouse 1 shock_post/shock_pre safe_post/safe_pre
% mouse 2 shock_post/shock_pre safe_post/safe_pre
% ...

ratio_fentry = nan(length(Dir.path), 2, 2);
diff_fentry = nan(length(Dir.path), 2, 2);

for i=1:length(Dir.path)  % fix mouse
    for j = 1:2 % fix either Pre/Post val, or Shock/Safe val
        for k = 1:2 % fix column of Pre/Post in 1st val or Shock/Safe in 2nd val
            if j == 1
                ratio_fentry(i, k, j) = EntryTime_mean(i, 1, k)/EntryTime_mean(i, 2, k); 
                diff_fentry(i, k, j) = EntryTime_mean(i, 1, k) - EntryTime_mean(i, 2, k);                
            elseif j == 2
                ratio_fentry(i, k, j) = EntryTime_mean(i, k, 2)/EntryTime_mean(i, k, 1);
                diff_fentry(i, k, j) = EntryTime_mean(i, k, 2) - EntryTime_mean(i, k, 1);
            end
        end
    end
end

%% Calculate number of entries into the shock/safe zone

NumEntries_shock = nan(length(Dir.path), num_sess, 2); % 4 tests, Pre and Post
NumEntries_safe = nan(length(Dir.path), num_sess, 2); % 4 tests, Pre and Post

for i = 1:length(a)
    for k=1:Test_iterr %iterator for tests
        temp = CalculateNumEntriesZone(a{i}.behavResources(id_Pre{i}(k))); % Each iterration - a certain trial for each zone
        NumEntries_shock(i,k,1) = temp(1); % shock zone pretest
        NumEntries_safe(i,k,1) = temp(2); % safe zone pretest
        temp = CalculateNumEntriesZone(a{i}.behavResources(id_Post{i}(k)));
        NumEntries_shock(i,k,2) = temp(1);% shock zone post
        NumEntries_safe(i,k,2) = temp(2);% safe zone post
    end  
end

NumEntries_mean = nan(length(Dir.path), 2, 2);
NumEntries_std = nan(length(Dir.path), 2, 2);
for izone = 1:2 % 1st val codes preTest, 2nd val codes postTest
    NumEntries_mean(:, 1, izone) = mean(squeeze(NumEntries_shock(:,:,izone)),2);
    NumEntries_mean(:, 2, izone) = mean(squeeze(NumEntries_safe(:,:,izone)),2);
    
    NumEntries_std(:, 1, izone) = std(squeeze(NumEntries_shock(:,:,izone)),0,2);
    NumEntries_std(:, 2, izone) = std(squeeze(NumEntries_safe(:,:,izone)),0,2);
end

%% Calculate ratios of number of entries into the shock/safe zone
% structure of ratio_nentries:
% val(:, :, 1) = 
%             Pre       Post
% mouse 1 shock/safe shock/safe
% mouse 2 shock/safe shock/safe
% ...
% 
% val(:, :, 2) = 
%                 Shock               Safe
% mouse 1 shock_post/shock_pre safe_post/safe_pre
% mouse 2 shock_post/shock_pre safe_post/safe_pre
% ...

ratio_nentries = nan(length(Dir.path), 2, 2);
diff_nentries = nan(length(Dir.path), 2, 2);
for i=1:length(Dir.path)  % fix mouse
    for j = 1:2 % fix either Pre/Post val, or Shock/Safe val
        for k = 1:2 % fix column of Pre/Post in 1st val or Shock/Safe in 2nd val
            if j == 1
                ratio_nentries(i, k, j) = NumEntries_mean(i, 1, k)/NumEntries_mean(i, 2, k); 
                diff_nentries(i, k, j) = NumEntries_mean(i, 1, k) - NumEntries_mean(i, 2, k);                 
            elseif j == 2
                ratio_nentries(i, k, j) = NumEntries_mean(i, k, 2)/NumEntries_mean(i, k, 1);
                diff_nentries(i, k, j) = NumEntries_mean(i, k, 2) - NumEntries_mean(i, k, 1);
            end
        end
    end
end

%% Calculate average speed in the shock and safe zone
% I skip the last point in ZoneIndices because length(Xtsd)=length(Vtsd)+1
% - UPD 18/07/2018 - Could do length(Start(ZoneEpoch))

Speed_shock = nan(length(Dir.path), num_sess, 2); % 4 tests Pre and Post
Speed_safe = nan(length(Dir.path), num_sess, 2); % 4 tests Pre and Post

for i = 1:length(a)
    for k=1:Test_iterr % iterator for tests
        temp = CalculateSpeedZone(a{i}.behavResources(id_Pre{i}(k)));
        Speed_shock(i,k,1) = temp(1);
        Speed_safe(i,k,1) = temp(2);
        
        temp = CalculateSpeedZone(a{i}.behavResources(id_Post{i}(k)));
        Speed_shock(i,k,2) = temp(1);
        Speed_safe(i,k,2) = temp(2);
    end 
end

Speed_mean = nan(length(Dir.path), 2, 2);
Speed_std = nan(length(Dir.path), 2, 2);

for izone = 1:2 % 1st val codes preTest, 2nd val codes postTest
    Speed_mean(:, 1, izone) = mean(squeeze(Speed_shock(:,:,izone)),2);
    Speed_mean(:, 2, izone) = mean(squeeze(Speed_safe(:,:,izone)),2);
    
    Speed_std(:, 1, izone) = std(squeeze(Speed_shock(:,:,izone)), 0, 2);
    Speed_std(:, 2, izone) = std(squeeze(Speed_safe(:,:,izone)), 0, 2);
end

%% Calculate speed distribution
% for izone = 1:2
%     for imouse = 1:length(a)
%         i = 1;
%         if izone == 1
%             for id_test = id_Pre{imouse}
%                 Speed_data(imouse).data_PreTest{i} = Data(a{imouse}.behavResources(id_test).CleanVtsd);
%                 i = i + 1;
%             end
%         else
%             for id_test = id_Post{imouse}
%                 Speed_data(imouse).data_PostTest{i} = Data(a{imouse}.behavResources(id_test).CleanVtsd);
%                 i = i + 1;
%             end
%         end
%     end
% end
% 
% conc_Pre_Speed = [Speed_data.data_PreTest];
% conc_Post_Speed = [Speed_data.data_PostTest];
% 
% for iconc = 1:length(conc_Pre_Speed)
%     length_conc_pre(iconc) = [length(conc_Pre_Speed{iconc})];
%     length_conc_post(iconc) = [length(conc_Post_Speed{iconc})];
% end
% length_conc_pre = length_conc_pre';
% max_length_conc_pre = max(length_conc_pre);
% length_conc_post = length_conc_post';
% max_length_conc_post = max(length_conc_post);
% 
% for iconc = 1:length(conc_Pre_Speed)
%     if length(conc_Pre_Speed{iconc}) == max_length_conc_pre
%         i = i + 1;
%     else
%         final_pre_array{iconc} = padarray(conc_Pre_Speed{iconc}, max_length_conc_pre - length(conc_Pre_Speed{iconc}), NaN, 'post');
%     end
%     if length(conc_Post_Speed{iconc}) == max_length_conc_post
%         i = i + 1;
%     else
%         final_post_array{iconc} = padarray(conc_Post_Speed{iconc}, max_length_conc_post - length(conc_Post_Speed{iconc}), NaN, 'post');
%     end
% end
% 
% final_final_post_array = [final_post_array{1} final_post_array{2}  final_post_array{3}  final_post_array{4} final_post_array{5} final_post_array{6} final_post_array{7} final_post_array{8} final_post_array{9} final_post_array{10} final_post_array{11} final_post_array{12} final_post_array{13} final_post_array{14} final_post_array{15} final_post_array{16} final_post_array{17} final_post_array{18} final_post_array{19} final_post_array{20} final_post_array{21} final_post_array{22} final_post_array{23} final_post_array{24} final_post_array{25} final_post_array{26} final_post_array{27} final_post_array{28}];
% final_final_pre_array = [final_pre_array{1} final_pre_array{2}  final_pre_array{3}  final_pre_array{4} final_pre_array{5} final_pre_array{6} final_pre_array{7} final_pre_array{8} final_pre_array{9} final_pre_array{10} final_pre_array{11} final_pre_array{12} final_pre_array{13} final_pre_array{14} final_pre_array{15} final_pre_array{16} final_pre_array{17} final_pre_array{18} final_pre_array{19} final_pre_array{20} final_pre_array{21} final_pre_array{22} final_pre_array{23} final_pre_array{24} final_pre_array{25} final_pre_array{26} final_pre_array{27} final_pre_array{28}];
% 
% final_final_post_array = final_post_array{1};
% for i=2:length(final_final_post_array)
%     final_final_post_array = [final_final_post_array; final_post_array{i}];
% end

% figure
% % ylim_out = SelectYlim(fhandle);
% subplot(121)
% histogram(final_final_pre_array);
% xlim([0, 20])
% ylim([0, 2.5e4])
% % ylim([ylim_out])
% xlabel('Speed (cm/s)')
% title('PreTests')
% makepretty
% 
% subplot(122)
% histogram(final_final_post_array);
% xlim([0, 20])
% ylim([0, 2.5e4])
% % ylim([ylim_out])
% xlabel('Speed (cm/s)')
% title('PostTests')
% makepretty

%% Calculate ratios of speed

% structure of ratio_speed:
% val(:, :, 1) = 
%             Pre       Post
% mouse 1 shock/safe shock/safe
% mouse 2 shock/safe shock/safe
% ...
% 
% val(:, :, 2) = 
%                 Shock               Safe
% mouse 1 shock_post/shock_pre safe_post/safe_pre
% mouse 2 shock_post/shock_pre safe_post/safe_pre
% ...

ratio_speed = nan(length(Dir.path), 2, 2);
diff_speed = nan(length(Dir.path), 2, 2);

for i=1:length(Dir.path)  % fix mouse
    for j = 1:2 % fix either Pre/Post val, or Shock/Safe val
        for k = 1:2 % fix column of Pre/Post in 1st val or Shock/Safe in 2nd val
            if j == 1
                ratio_speed(i, k, j) = Speed_mean(i, 1, k)/Speed_mean(i, 2, k); 
                diff_speed(i, k, j) = Speed_mean(i, 1, k) - Speed_mean(i, 2, k); 
            elseif j == 2
                ratio_speed(i, k, j) = Speed_mean(i, k, 2)/Speed_mean(i, k, 1);
                diff_speed(i, k, j) = Speed_mean(i, k, 2) - Speed_mean(i, k, 1);
            end
        end
    end
end

%% Plot shock zone
% Axes
fh = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);
Occupancy_Axes = axes('position', [0.07 0.55 0.41 0.41]);
NumEntr_Axes = axes('position', [0.55 0.55 0.41 0.41]);
First_Axes = axes('position', [0.07 0.05 0.41 0.41]);
Speed_Axes = axes('position', [0.55 0.05 0.41 0.41]);

% Occupancy
axes(Occupancy_Axes);
[p_occ,h_occ, her_occ] = PlotErrorBarN_DB([occup_mean(:,1,1)*100 occup_mean(:,1,2)*100],...
    'barcolors', [.8 .9 .7], 'barwidth', 0.6, 'newfig', 0, 'showpoints', SP);
h_occ.FaceColor = 'flat';
h_occ.CData(2,:) = [0 .5 .3];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_occ, 'LineWidth', 3);
set(her_occ, 'LineWidth', 3);
line(xlim,[21.5 21.5],'Color','k','LineStyle','--','LineWidth',3);
ylabel('% time');
title('Occupancy percentage of shock zone', 'FontSize', 14);
ylim([0 100])

axes(NumEntr_Axes);
[p_nent,h_nent, her_nent] = PlotErrorBarN_DB([NumEntries_mean(:,1,1) NumEntries_mean(:,1,2)],...
    'barcolors', [.7 .7 .4], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_nent.FaceColor = 'flat';
h_nent.CData(2,:) = [.8 .7 .3];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_nent, 'LineWidth', 3);
set(her_nent, 'LineWidth', 3);
ylabel('Number of entries');
title('Number of entries to the shock zone', 'FontSize', 14);
ylim([0 9])

axes(First_Axes);
[p_first,h_first, her_first] = PlotErrorBarN_DB([EntryTime_mean(:,1,1) EntryTime_mean(:,1,2)],...
    'barcolors', [.2 .6 .7], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_first.FaceColor = 'flat';
h_first.CData(2,:) = [0 .4 .5];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 18, 'FontWeight', 'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_first, 'LineWidth', 3);
set(her_first, 'LineWidth', 3);
ylabel('Time (s)');
title('First time to enter the shock zone', 'FontSize', 14);
ylim([0 250])

axes(Speed_Axes);
[p_speed,h_speed, her_speed] = PlotErrorBarN_DB([Speed_mean(:,1,1) Speed_mean(:,1,2)],...
    'barcolors', [.7 .3 .3], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_speed.FaceColor = 'flat';
h_speed.CData(2,:) = [.7 0 .1];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_speed, 'LineWidth', 3);
set(her_speed, 'LineWidth', 3);
ylabel('Speed (cm/s)');
title('Average speed in the shock zone', 'FontSize', 14);
ylim([0 6])

%% Plot safe zone
% Axes
fh = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);
Occupancy_Axes = axes('position', [0.07 0.55 0.41 0.41]);
NumEntr_Axes = axes('position', [0.55 0.55 0.41 0.41]);
First_Axes = axes('position', [0.07 0.05 0.41 0.41]);
Speed_Axes = axes('position', [0.55 0.05 0.41 0.41]);

% Occupancy
axes(Occupancy_Axes);
[p_occ,h_occ, her_occ] = PlotErrorBarN_DB([occup_mean(:,2,1)*100 occup_mean(:,2,2)*100],...
    'barcolors', [.8 .9 .7], 'barwidth', 0.6, 'newfig', 0, 'showpoints', SP);
h_occ.FaceColor = 'flat';
h_occ.CData(2,:) = [0 .5 .3];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_occ, 'LineWidth', 3);
set(her_occ, 'LineWidth', 3);
line(xlim,[21.5 21.5],'Color','k','LineStyle','--','LineWidth',3);
ylabel('% time');
title('Occupancy percentage of safe zone', 'FontSize', 14);
ylim([0 100])

axes(NumEntr_Axes);
[p_nent,h_nent, her_nent] = PlotErrorBarN_DB([NumEntries_mean(:,2,1) NumEntries_mean(:,2,2)],...
    'barcolors', [.7 .7 .4], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_nent.FaceColor = 'flat';
h_nent.CData(2,:) = [.8 .7 .3];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_nent, 'LineWidth', 3);
set(her_nent, 'LineWidth', 3);
ylabel('Number of entries');
title('Number of entries to the safe zone', 'FontSize', 14);
ylim([0 9])

axes(First_Axes);
[p_first,h_first, her_first] = PlotErrorBarN_DB([EntryTime_mean(:,2,1) EntryTime_mean(:,2,2)],...
    'barcolors', [.2 .6 .7], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_first.FaceColor = 'flat';
h_first.CData(2,:) = [0 .4 .5];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 18, 'FontWeight', 'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_first, 'LineWidth', 3);
set(her_first, 'LineWidth', 3);
ylabel('Time (s)');
title('First time to enter the safe zone', 'FontSize', 14);
ylim([0 250])

axes(Speed_Axes);
[p_speed,h_speed, her_speed] = PlotErrorBarN_DB([Speed_mean(:,2,1) Speed_mean(:,2,2)],...
    'barcolors', [.7 .3 .3], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_speed.FaceColor = 'flat';
h_speed.CData(2,:) = [.7 0 .1];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_speed, 'LineWidth', 3);
set(her_speed, 'LineWidth', 3);
ylabel('Speed (cm/s)');
title('Average speed in the safe zone', 'FontSize', 14);
ylim([0 6])

%% Plot occupancy and speed ratios
% % Axes
% fr_occup_speed = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);
% 
% Ratio_occup_prepost_Axes = axes('position', [0.08 0.55 0.35 0.35]);
% Ratio_occup_shocksafe_Axes = axes('position', [0.58 0.55 0.35 0.35]);
% 
% Ratio_speed_prepost_Axes = axes('position', [0.08 0.05 0.35 0.35]);
% Ratio_speed_shocksafe_Axes = axes('position', [0.58 0.05 0.35 0.35]);
% 
% % Occupancy
% axes(Ratio_occup_prepost_Axes);
% [p_occ_pp,h_occ_pp, her_occ_pp] = PlotErrorBarN_DB([ratio_occup(:,1,1) ratio_occup(:,2,1)],...
%     'barcolors', [1 .8 .8], 'barwidth', 0.4, 'newfig', 0, 'showpoints', SP);
% h_occ_pp.FaceColor = 'flat';
% h_occ_pp.CData(2,:) = [0 .5 .5];
% set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
% set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(gca, 'LineWidth', 2);           
% set(h_occ_pp, 'LineWidth', 2);
% set(her_occ_pp, 'LineWidth', 2);
% line(xlim,[1 1],'Color','k','LineStyle','--','LineWidth',1);
% ylabel('Shock/Safe');
% title('Shock/Safe ratio of Pre/Post occupancy', 'FontSize', 14);
% 
% axes(Ratio_occup_shocksafe_Axes);
% [p_occ_ss,h_occ_ss, her_occ_ss] = PlotErrorBarN_DB([ratio_occup(:,1,2) ratio_occup(:,2,2)],...
%     'barcolors', [1 .8 .8], 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
% h_occ_ss.FaceColor = 'flat';
% h_occ_ss.CData(2,:) = [0 .5 .5];
% set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
% set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(gca, 'LineWidth', 2);
% set(h_occ_ss, 'LineWidth', 2);
% set(her_occ_ss, 'LineWidth', 2);
% line(xlim,[1 1],'Color','k','LineStyle','--','LineWidth',1);
% ylabel('Post/Pre');
% title('Post/Pre ratio of shock/safe zones occupancy', 'FontSize', 14);
% 
% % Speed
% axes(Ratio_speed_prepost_Axes);
% [p_speed_pp,h_speed_pp, her_speed_pp] = PlotErrorBarN_DB([ratio_speed(:,1,1) ratio_speed(:,2,1)],...
%     'barcolors', [1 .8 .8], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
% h_speed_pp.FaceColor = 'flat';
% h_speed_pp.CData(2,:) = [0 .5 .5];
% set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
% set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(gca, 'LineWidth', 2);
% set(h_speed_pp, 'LineWidth', 2);
% set(her_speed_pp, 'LineWidth', 2);
% line(xlim,[1 1],'Color','k','LineStyle','--','LineWidth',1);
% ylabel('Shock/Safe');
% title('Shock/Safe ratio of Pre/Post speed', 'FontSize', 14);
% 
% axes(Ratio_speed_shocksafe_Axes);
% [p_speed_ss,h_speed_ss, her_speed_ss] = PlotErrorBarN_DB([ratio_speed(:,1,2) ratio_speed(:,2,2)],...
%     'barcolors', [1 .8 .8], 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
% h_speed_ss.FaceColor = 'flat';
% h_speed_ss.CData(2,:) = [0 .5 .5];
% set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
% set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(gca, 'LineWidth', 2);
% set(h_speed_ss, 'LineWidth', 2);
% set(her_speed_ss, 'LineWidth', 2);
% line(xlim,[1 1],'Color','k','LineStyle','--','LineWidth',1);
% ylabel('Post/Pre');
% title('Post/Pre ratio of shock/safe zones speed', 'FontSize', 14);

%% Plot occupancy and speed differences
% Axes
fd_occup_speed = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);

Diff_occup_prepost_Axes = axes('position', [0.08 0.55 0.35 0.35]);
Diff_occup_shocksafe_Axes = axes('position', [0.58 0.55 0.35 0.35]);

Diff_speed_prepost_Axes = axes('position', [0.08 0.05 0.35 0.35]);
Diff_speed_shocksafe_Axes = axes('position', [0.58 0.05 0.35 0.35]);

% Occupancy
axes(Diff_occup_prepost_Axes);
[p_occ_pp,h_occ_pp, her_occ_pp] = PlotErrorBarN_DB([diff_occup(:,1,1) diff_occup(:,2,1)],...
    'barcolors', [.8 .9 .7], 'barwidth', 0.4, 'newfig', 0, 'showpoints', SP);
h_occ_pp.FaceColor = 'flat';
h_occ_pp.CData(2,:) = [0 .5 .3];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);           
set(h_occ_pp, 'LineWidth', 2);
set(her_occ_pp, 'LineWidth', 2);
ylabel('Shock - Safe');
title('Shock - Safe difference of Pre/Post occupancy', 'FontSize', 14);
t1 = text(.8, .75, ['Shock > Safe'], 'sc');
t2 = text(.8, .25, ['Shock < Safe'], 'sc');

axes(Diff_occup_shocksafe_Axes);
[p_occ_ss,h_occ_ss, her_occ_ss] = PlotErrorBarN_DB([diff_occup(:,1,2) diff_occup(:,2,2)],...
    'barcolors', [.8 .9 .7], 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_occ_ss.FaceColor = 'flat';
h_occ_ss.CData(2,:) = [0 .5 .3];
set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);
set(h_occ_ss, 'LineWidth', 2);
set(her_occ_ss, 'LineWidth', 2);
ylabel('Post - Pre');
title('Post - Pre difference of shock/safe zones occupancy', 'FontSize', 14);
t1 = text(.8, .8, ['Post > Pre'], 'sc');
t2 = text(.8, .3, ['Post < Pre'], 'sc');

% Speed
axes(Diff_speed_prepost_Axes);
[p_speed_pp,h_speed_pp, her_speed_pp] = PlotErrorBarN_DB([diff_speed(:,1,1) diff_speed(:,2,1)],...
    'barcolors', [.7 .3 .3], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_speed_pp.FaceColor = 'flat';
h_speed_pp.CData(2,:) = [.7 0 .1];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);
set(h_speed_pp, 'LineWidth', 2);
set(her_speed_pp, 'LineWidth', 2);
ylabel('Shock - Safe');
title('Shock - Safe difference of Pre/Post speed', 'FontSize', 14);
t1 = text(.8, .75, ['Shock > Safe'], 'sc');
t2 = text(.8, .25, ['Shock < Safe'], 'sc');

axes(Diff_speed_shocksafe_Axes);
[p_speed_ss,h_speed_ss, her_speed_ss] = PlotErrorBarN_DB([diff_speed(:,1,2) diff_speed(:,2,2)],...
    'barcolors', [.7 .3 .3], 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_speed_ss.FaceColor = 'flat';
h_speed_ss.CData(2,:) = [.7 0 .1];
set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);
set(h_speed_ss, 'LineWidth', 2);
set(her_speed_ss, 'LineWidth', 2);
ylabel('Post - Pre');
title('Post - Pre difference of shock/safe zones speed', 'FontSize', 14);
t1 = text(.8, .875, ['Post > Pre'], 'sc');
t2 = text(.8, .3, ['Post < Pre'], 'sc');

%% Plot entries ratios
% % Axes
% fr_entries = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);
% 
% Ratio_nentries_prepost_Axes = axes('position', [0.08 0.55 0.35 0.35]);
% Ratio_nentries_shocksafe_Axes = axes('position', [0.58 0.55 0.35 0.35]);
% 
% Ratio_fentry_prepost_Axes = axes('position', [0.08 0.05 0.35 0.35]);
% Ratio_fentry_shocksafe_Axes = axes('position', [0.58 0.05 0.35 0.35]);
% 
% % Number of entries
% axes(Ratio_nentries_prepost_Axes);
% [p_nentries_pp,h_nentries_pp, her_nentries_pp] = PlotErrorBarN_DB([ratio_nentries(:,1,1) ratio_nentries(:,2,1)],...
%     'barcolors', [1 .8 .8], 'barwidth', 0.4, 'newfig', 0, 'showpoints', SP);
% h_nentries_pp.FaceColor = 'flat';
% h_nentries_pp.CData(2,:) = [0 .5 .5];
% set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
% set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(gca, 'LineWidth', 2);           
% set(h_nentries_pp, 'LineWidth', 2);
% set(her_nentries_pp, 'LineWidth', 2);
% line(xlim,[1 1],'Color','k','LineStyle','--','LineWidth',1);
% ylabel('Shock/Safe');
% title('Shock/Safe ratio of Pre/Post number of entries to each zone', 'FontSize', 14);
% 
% axes(Ratio_nentries_shocksafe_Axes);
% [p_nentries_ss,h_nentries_ss, her_nentries_ss] = PlotErrorBarN_DB([ratio_nentries(:,1,2) ratio_nentries(:,2,2)],...
%     'barcolors', [1 .8 .8], 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
% h_nentries_ss.FaceColor = 'flat';
% h_nentries_ss.CData(2,:) = [0 .5 .5];
% set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
% set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(gca, 'LineWidth', 2);
% set(h_nentries_ss, 'LineWidth', 2);
% set(her_nentries_ss, 'LineWidth', 2);
% line(xlim,[1 1],'Color','k','LineStyle','--','LineWidth',1);
% ylabel('Post/Pre');
% title('Post/Pre ratio of shock/safe zones entries', 'FontSize', 14);
% 
% % First entry ratios
% axes(Ratio_fentry_prepost_Axes);
% [p_fentry_pp,h_fentry_pp, her_fentry_pp] = PlotErrorBarN_DB([ratio_fentry(:,1,1) ratio_fentry(:,2,1)],...
%     'barcolors', [1 .8 .8], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
% h_fentry_pp.FaceColor = 'flat';
% h_fentry_pp.CData(2,:) = [0 .5 .5];
% set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
% set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(gca, 'LineWidth', 2);
% set(h_fentry_pp, 'LineWidth', 2);
% set(her_fentry_pp, 'LineWidth', 2);
% line(xlim,[1 1],'Color','k','LineStyle','--','LineWidth',1);
% ylabel('Shock/Safe');
% title('Shock/Safe ratio of the time of the first entry to the zone during Pre/Post', 'FontSize', 14);
% 
% axes(Ratio_fentry_shocksafe_Axes);
% [p_fentry_ss,h_fentry_ss, her_fentry_ss] = PlotErrorBarN_DB([ratio_fentry(:,1,2) ratio_fentry(:,2,2)],...
%     'barcolors', [1 .8 .8], 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
% h_fentry_ss.FaceColor = 'flat';
% h_fentry_ss.CData(2,:) = [0 .5 .5];
% set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
% set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(gca, 'LineWidth', 2);
% set(h_fentry_ss, 'LineWidth', 2);
% set(her_fentry_ss, 'LineWidth', 2);
% line(xlim,[1 1],'Color','k','LineStyle','--','LineWidth',1);
% ylabel('Post/Pre');
% title('Post/Pre ratio of first entries to shock/safe zones', 'FontSize', 14);

%% Plot entries differences
% Axes
fd_entries = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);

Diff_nentries_prepost_Axes = axes('position', [0.08 0.55 0.35 0.35]);
Diff_nentries_shocksafe_Axes = axes('position', [0.58 0.55 0.35 0.35]);

Diff_fentry_prepost_Axes = axes('position', [0.08 0.05 0.35 0.35]);
Diff_fentry_shocksafe_Axes = axes('position', [0.58 0.05 0.35 0.35]);

% Number of entries
axes(Diff_nentries_prepost_Axes);
[p_nentries_pp,h_nentries_pp, her_nentries_pp] = PlotErrorBarN_DB([diff_nentries(:,1,1) diff_nentries(:,2,1)],...
    'barcolors', [.7 .7 .4], 'barwidth', 0.4, 'newfig', 0, 'showpoints', SP);
h_nentries_pp.FaceColor = 'flat';
h_nentries_pp.CData(2,:) = [.8 .7 .3];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);           
set(h_nentries_pp, 'LineWidth', 2);
set(her_nentries_pp, 'LineWidth', 2);
ylabel('Shock - Safe');
title('Shock - Safe difference of Pre/Post number of entries to each zone', 'FontSize', 14);
t1 = text(.8, .9, ['Shock > Safe'], 'sc');
t2 = text(.8, .3, ['Shock < Safe'], 'sc');
ylim([-3, 3]);

axes(Diff_nentries_shocksafe_Axes);
[p_nentries_ss,h_nentries_ss, her_nentries_ss] = PlotErrorBarN_DB([diff_nentries(:,1,2) diff_nentries(:,2,2)],...
    'barcolors', [.7 .7 .4], 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_nentries_ss.FaceColor = 'flat';
h_nentries_ss.CData(2,:) = [.8 .7 .3];
set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);
set(h_nentries_ss, 'LineWidth', 2);
set(her_nentries_ss, 'LineWidth', 2);
ylabel('Post - Pre');
title('Post - Pre difference of shock/safe zones entries', 'FontSize', 14);
t1 = text(.8, .85, ['Post > Pre'], 'sc');
t2 = text(.8, .35, ['Post < Pre'], 'sc');

% First entry differences
axes(Diff_fentry_prepost_Axes);
[p_fentry_pp,h_fentry_pp, her_fentry_pp] = PlotErrorBarN_DB([diff_fentry(:,1,1) diff_fentry(:,2,1)],...
    'barcolors', [.2 .6 .7], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_fentry_pp.FaceColor = 'flat';
h_fentry_pp.CData(2,:) = [0 .4 .5];
set(gca,'Xtick',[1:2],'XtickLabel',{'Pre', 'Post'});
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);
set(h_fentry_pp, 'LineWidth', 2);
set(her_fentry_pp, 'LineWidth', 2);
ylabel('Shock - Safe');
title('Shock - Safe difference of the time of the first entry to the zone during Pre/Post', 'FontSize', 14);
t1 = text(.8, .65, ['Shock > Safe'], 'sc');
t2 = text(.8, .15, ['Shock < Safe'], 'sc');

axes(Diff_fentry_shocksafe_Axes);
[p_fentry_ss,h_fentry_ss, her_fentry_ss] = PlotErrorBarN_DB([diff_fentry(:,1,2) diff_fentry(:,2,2)],...
    'barcolors', [.2 .6 .7], 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_fentry_ss.FaceColor = 'flat';
h_fentry_ss.CData(2,:) = [0 .4 .5];
set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);
set(h_fentry_ss, 'LineWidth', 2);
set(her_fentry_ss, 'LineWidth', 2);
ylabel('Post - Pre');
title('Post - Pre difference of first entries to shock/safe zones', 'FontSize', 14);
t1 = text(.8, .6, ['Post > Pre'], 'sc');
t2 = text(.8, .1, ['Post < Pre'], 'sc');

%% Save
if IsSave
    saveas(fh,[pwd '/' FigName '.fig']);
    saveFigure(fh, pwd, dirsave);
    saveFigure(fr_occup_speed, pwd, dirsave);
    saveFigure(fr_entries, pwd, dirsave);
    saveFigure(fd_occup_speed, pwd, dirsave);
    saveFigure(fd_entries, pwd, dirsave);
end

%% Save ratios

if save_res == 1
    save([foldertosave filesep 'Behav_ratios.mat'], 'ratio_occup', 'diff_occup', 'ratio_fentry', 'diff_fentry', 'ratio_nentries', 'diff_nentries', 'ratio_speed', 'diff_speed');
end

%% Testing adequacy of nentries
% Test Pre trajectories
% for i = 3:6
%     figure, plot(Data(a{1}.behavResources(i).CleanAlignedXtsd), Data(a{1}.behavResources(i).CleanAlignedYtsd))
% end
% % % Test Post trajectories
% for i = 12:15
%     figure, plot(Data(a{1}.behavResources(i).CleanAlignedXtsd), Data(a{1}.behavResources(i).CleanAlignedYtsd))
% end

end