function id_TT = PyramLayerSorting(RippleGroups, TT) 
%%
% This function filters neurons from pyramidal layer of hippocampus with
% certain indixes which are indicated in SpikeData.mat -> 'RippleGroups'.
%%
id_TT = zeros(length(TT), 1); % Preallocation of the array
for iRG = 1:length(RippleGroups)
    for iTT = 1:length(TT)
        if TT{iTT}(1) == RippleGroups(iRG)
            id_TT(iTT) = iTT;
        end
    end
end
id_TT = nonzeros(id_TT);

end
