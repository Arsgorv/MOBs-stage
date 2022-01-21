AddMyPaths_Arsenii

% Changable parameters
nmice = [1199];
experiment = "PAG";
TemplatePeriod = {'condRip'};% 'wake' 'cond', 'condFree', 'postRip', 'condRip'
numtemplates = 13; %Number of PCs.
savefig = 1; %if you want to save figures to your pathway - put 1
binsize = 0.025*1e4; %0.1*1e4 = 100 ms; 0.025*1e4 = 25 ms.
calc = "PCA"; % TRUE if you want to calculate PCA (pcacov)
issavefig = 0; 

% close all;

for i=1:length(TemplatePeriod)
    [PCAtemplates, EV_PCA] = react_pca_ica_AG(nmice,experiment,TemplatePeriod{1,i}, numtemplates, "PCA", 'binsize', binsize, 'issavefig', issavefig);

    EV_PCA = [EV_PCA{:}];
    EV = sum(EV_PCA(1:numtemplates))

    PCAtemplates = [PCAtemplates{:}];
    PCAtemplates = PCAtemplates(:,1:numtemplates);

    ICAtemplates = react_pca_ica_AG(nmice,experiment,TemplatePeriod{1,i}, numtemplates, "ICA", 'binsize', binsize, 'issavefig', issavefig);
    ICAtemplates = [ICAtemplates{:}];

    
    corrcoef(PCAtemplates(:,1), ICAtemplates(:,1));

    PC_correl = zeros(numtemplates);

    for nPC_PCA = 1:numtemplates
        for nPC_ICA = 1:numtemplates
            A = corrcoef(PCAtemplates(:,nPC_PCA), ICAtemplates(:, nPC_ICA));
            PC_correl(nPC_PCA, nPC_ICA) = A(1,2);
        end
    end
end

row = strings([1, numtemplates]);
col = strings([1, numtemplates]);
for nPC_PCA = 1:numtemplates
    row(nPC_PCA) = strjoin({'PCA' num2str(nPC_PCA)});
    col(nPC_PCA) = strjoin({'ICA' num2str(nPC_PCA)});
end

imagesc(PC_correl)
set(gca,'Xtick',1:numtemplates,'XTickLabel',col);
set(gca, 'Ytick', 1:numtemplates, 'YtickLabel', row);
caxis([-1 1]);
colormap jet;
colorbar;
title("Correlation between ICA and PCA components");

if issavefig == 1
     saveas(gcf, ['PCA_ICA_Correlations' nmice(1) '_templateEpoch_' TemplatePeriod{1} '.png']);
end

% close all;

num_sign_PCA = [];
neuron_sign_PCA = {};
temp_neuron_PCA = [];
num_sign_ICA = [];
neuron_sign_ICA = {};
temp_neuron_ICA = [];

for temp=1:numtemplates
    temporary_PCA = 0;
    temporary_ICA = 0;
    for neuron = 1:size(PCAtemplates,1)
        if abs(PCAtemplates(neuron,temp)) > (mean(PCAtemplates(:,temp)) + 2*std(PCAtemplates(:,temp)))
            temporary_PCA = temporary_PCA + 1;
            temp_neuron_PCA = [temp_neuron_PCA, neuron];
        end

        if abs(ICAtemplates(neuron,temp)) > (mean(ICAtemplates(:,temp)) + 2*std(ICAtemplates(:,temp)))
            temporary_ICA = temporary_ICA + 1;
            temp_neuron_ICA = [temp_neuron_ICA, neuron];
        end

    end 
    num_sign_PCA(temp) = temporary_PCA;
    neuron_sign_PCA{temp} = [temp_neuron_PCA];
    num_sign_ICA(temp) = temporary_ICA;
    neuron_sign_ICA{temp} = [temp_neuron_ICA];
end

p = ranksum(num_sign_PCA(:), num_sign_ICA(:));

figure
boxplot([num_sign_PCA(:), num_sign_ICA(:)], 'Labels', {'PCA', 'ICA'}, 'Whisker',1)
text(1, 1, num2str(p));
title('Comparison of mean number of significant neurons per cell assemblies')

if issavefig == 1
    saveas(gcf, ['PCA_ICA_CompSignW_param_' nmice(1) '_templateEpoch_' TemplatePeriod{1} '.png']);
end

% close all;
