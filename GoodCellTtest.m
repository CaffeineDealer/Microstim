function [p] = GoodCellTtest(frst,frbl,ch)


stim = squeeze(frst(ch,:,:,:));
baseline = squeeze(frbl(ch,:,:,:));
[h,p] = ttest(baseline(:),stim(:));
if p <= 0.05
    disp('Yo! You got a Cell!')
elseif p > 0.05
    warning('Not a Cell, Life Sucks!')
end