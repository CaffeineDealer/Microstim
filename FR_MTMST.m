function [FRmean FRmst err] = FR_MTMST(datapath,opt,type,flag,MS)





inf = load([datapath sprintf('N%sMS.mat',type)]); 
MST = inf.info; 
clear info inf
inf = load([datapath sprintf('%s%s.mat',opt,type)]); 
MT = inf.info; 
clear info inf
if MS == 1
    MT(7) = [];
    MT(4) = []; 
    MT(1) = [];
end

c = 1;
for i = 1:size(MT,2)
    fr = MT(i).FR(:,1) ./ max(MT(i).FR(:,1));
    err(i,1) = std(fr) / sqrt(size(fr,1));
    FRmean(i,1) = mean(fr);
    z = size(MST(i).FR(:,1),1);
    FRmst(c:c+z-1,1) = MST(i).FR(:,1);
    FRmst(c:c+z-1,2) = i;
    c = c + z;
end

FRmst(:,1) = FRmst(:,1) ./ max(FRmst(:,1));

