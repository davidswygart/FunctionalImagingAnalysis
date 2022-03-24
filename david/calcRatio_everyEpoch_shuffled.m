function [dF_SI, dFoF_SI, df_perm, dfof_perm] = calcRatio_everyEpoch_shuffled(smallG, bigG, smallT, bigT)
%% average accross time for each epoch
[smallGPre, smallGStim] = pullPreStim(smallG,smallT);
smallGPre = squeeze(mean(smallGPre, 3, 'omitnan'));
smallGStim = squeeze(mean(smallGStim, 3, 'omitnan'));

[bigGPre, bigGStim] = pullPreStim(bigG,bigT);
bigGPre = squeeze(mean(bigGPre, 3, 'omitnan'));
bigGStim = squeeze(mean(bigGStim, 3, 'omitnan'));

%% calculate responses (dF and dFoF) for each epoch
dF_small = smallGStim - smallGPre;
dF_big = bigGStim - bigGPre;

dFoF_small = dF_small ./ smallGPre;
dFoF_big = dF_big ./ bigGPre;

%% calculate SI
dF_SI = 1 - dF_big ./ dF_small;
dFoF_SI = 1 - dFoF_big ./ dFoF_small;

%% average SI accross epochs
dF_SI = median(dF_SI,3);
dFoF_SI = median(dFoF_SI,3);


%% calculate SI for every permuation of the epochs

[nr,nc,ne] = size(dF_big);

if ne > 7
    [~, epochPerms] = sort(rand(1000,ne),2); %just do 5k permutations if more than 8 epochs
else
    epochPerms = perms(1:ne); % 7 epochs will have 5k permutations
end

numPerms = size(epochPerms,1);

perm_sup_df = nan(nr,nc, ne*numPerms);
perm_sup_dfof = nan(nr,nc, ne*numPerms);

for i=1:numPerms
    inds = epochPerms(i,:);
    perm_sup_df(:,:,(i*ne-ne+1):(i*ne)) = 1 - dF_big ./ dF_small(:,:,inds);
    perm_sup_dfof(:,:,(i*ne-ne+1):(i*ne)) = 1 - dFoF_big ./ dFoF_small(:,:,inds);
end

 df_perm = median(perm_sup_df,3);
 dfof_perm = median(perm_sup_dfof,3);
end