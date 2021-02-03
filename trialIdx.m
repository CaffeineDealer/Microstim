function [comb] = trialIdx(info)





comb = [];
comb.allcomb = allcomb([unique(info(:,3))],[unique(info(:,4))],[unique(info(:,7))]);
for i = 1:size(comb.allcomb,1)
    comb.comb(i).idx = info(info(:,3)==comb.allcomb(i,1) & info(:,4)==comb.allcomb(i,2) & info(:,7)==comb.allcomb(i,3),12);
end