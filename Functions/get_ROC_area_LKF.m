 %%             get_ROC_area_LKF.m  
%%  function [AUC, ROC_rate] = get_ROC_area_LKF(spk_count1,spk_count2)
%%  modified from get_roc_curve.m
%%  input:  spike count in each trial, for each of two conditions
%       spk_count1  --- should be column vector, if you give row vector, the program will reshape it into column vector
%       spk_count2  --- should be column vector, if you give row vector, the program will reshape it into column vector
%                       usually spk_count1 is the distribution which has lower mean
%%   output:
%       AUC         ---  area under the curve
%       ROC_rate    --- first colum vector was false alarm rate, the second column was the hit rate
%% spk_count1 and spk_count2 should be column vector
%% spk_count1 is the distribution which has lower mean
%
%% for example, spike count in each trial 
% spk_count1=[     1
%      2
%      5
%      5
%      1
%      1
%      5
%      4
%      0
%      1];
% spk_count2 = [
%      6
%      6
%      5
%      1
%      2
%      8
%     10
%      5
%      2
%      1];
%      
function [AUC, ROC_rate] = get_ROC_area_LKF(spk_count1,spk_count2)
    if size(spk_count1,2)~=1&&size(spk_count1,1)==1
        spk_count1 = spk_count1';
    end
    
    if size(spk_count2,2)~=1&&size(spk_count2,1)==1
        spk_count2 = spk_count2';
    end
    
    n_thresh = 30;
    spk_index = [zeros(size(spk_count1));ones(size(spk_count2))];
    spk_c = [spk_count1;spk_count2];
    tmp = unique(spk_c);
    spk_range = linspace(tmp(1),tmp(end),n_thresh)';
    
    ROC_rate = zeros(n_thresh,2);
    for i = 1:n_thresh
        index = ge(spk_c,spk_range(i));
        
        ROC_rate(i,1) = sum(index==1&spk_index==0)/sum(spk_index==0);
        ROC_rate(i,2) = sum(index==1&spk_index==1)/sum(spk_index==1);
    end
    AUC = trapz(ROC_rate(end:-1:1,1),ROC_rate(end:-1:1,2));
    %ROC_rate = ROC_rate(:)';
    %output = [ROC_rate,AUC];
%     
%     figure;
% plot(ROC_rate(:,1), ROC_rate(:,2),'linewidth',3)
% x = xticks;
% hold on 
% plot(x,x,'k-')
end



%     Z = trapz(X,Y) computes the integral of Y with respect to X using
%     the trapezoidal method.  X and Y must be vectors of the same
%     length, or X must be a column vector and Y an array whose first
%     non-singleton dimension is length(X).