function [outputArg1] = calculate_dprime(nHit,nGo,nCR,nNogo)
%UNTITLED2 此处显示有关此函数的摘要
%   输入为hit/CR/Go/Nogo的数量
%   输出为dprime
if nHit == nGo
    hit_rate = 1-1/(2*nGo);
elseif nHit == 0
    hit_rate = 1/(2*nGo);
else
    hit_rate = nHit/nGo;
end
if nCR == nNogo
    FA_rate = 1/(2*nNogo);
elseif nCR == 0
    FA_rate = 1-1/(2*nNogo);
else
    FA_rate = 1-nCR/nNogo;
end
dprime = norminv(hit_rate)-norminv(FA_rate);

outputArg1 = dprime;
end

