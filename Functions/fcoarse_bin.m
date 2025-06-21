% function [cur_hist] = Fcoarse_bin(org_hist, bin_size)
% input: 
% org_hist - original histogram, bin at 0.1ms
% bin_size - the bin size of coarse bin
% output:
% cur_hist - the histogram bin at 'bin_size

function [cur_hist] = fcoarse_bin(org_hist, bin_size)

% coarse bin for 'org_hist'
tmp = size(org_hist);
tmp_1 = reshape(org_hist(:,1:bin_size*floor(tmp(2)/bin_size)),tmp(1),bin_size,floor(tmp(2)/bin_size));
cur_hist = sum(tmp_1,2);
cur_hist = squeeze(cur_hist);
%ave_cur_hist = mean(cur_hist,1);