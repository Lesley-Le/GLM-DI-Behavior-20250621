function [actualAUC,significant,significance] = Permutation_ROC_20250602(trial_type,signal)
%UNTITLED3 Calculate the AUC of 2 types of trials
%           and use permutation test to confirm whether this AUC is
%           significant
%   trial_type: include 2 trial types, must be marked by 1 & 2
%   signal: corresponding signals of the trial_type
%
%   Output: significant ---- 1 means the actual AUC is significant
%                       ---- 0 means not significant
%           significance ---- the exact accumulative probability of
%           actualAUC
%           actualAUC ---- the AUC calculated by the real trial type seq
%

shuffle_times = 1000;
alpha = 0.01;
Type = trial_type;
signal1 = signal(trial_type == 2);%the signal1 usually with lower mean
signal2 = signal(trial_type == 1);

realAUC = get_ROC_area_LKF(signal1,signal2);

%%Permutation Test: shuffle 1k times
for i = 1:shuffle_times
    index = randperm(length(Type));
    shuffle_data = Type(index);
    signal1 = signal(shuffle_data == 1);
    signal2 = signal(shuffle_data == 2);
    shuffle_DI(i,1) = get_ROC_area_LKF(signal1,signal2);
end

shuffle_DI_sorted = sort(shuffle_DI);

%%
delta = abs(shuffle_DI_sorted-realAUC);
[~,min_id] = min(delta);
significance = min_id/shuffle_times;
if realAUC<shuffle_DI_sorted(shuffle_times*alpha) || realAUC>=shuffle_DI_sorted(shuffle_times*(1-alpha))
    significant = 1;
else
    significant = 0;
end

actualAUC = realAUC;

end