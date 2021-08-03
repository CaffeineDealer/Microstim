function [rmean mmean] = StimRespMST(datapath,opt,type,MS)


inf = load([datapath sprintf('N%sMS.mat',type)]); 
MST = inf.info;
clear inf
inf = load([datapath sprintf('%s%s.mat',opt,type)]); 
MT = inf.info; 
clear inf
if MS == 1
    MT(7) = [];
    MT(4) = []; 
    MT(1) = [];
end


c = 1;
for i = 1:size(MST,2)
    for j = 1:size(MST(i).xr,2)
        mst = reshape(MST(i).xr(j).firing,72,1);
        mt = [];
        for z = 1:size(MT(i).xr,2)
            mt(:,z) = reshape(MT(i).xr(z).firing,72,1);
        end
        mtmst = [mst mt];
        mtmst = sortrows(mtmst,1);
        m(:,j) = mean(mtmst(:,2:end),2);
        r(:,j) = mtmst(:,1);
%         r(:,j) = sortrows(reshape(MST(i).xr(j).firing,72,1));
    end
    r = r ./ max(r);
    m = m ./ max(m);
    z = size(r,2);
    rall(:,c:c+z-1) = r;
    mall(:,c:c+z-1) = m;
    c = c + z;
    r = [];
    m = [];
    mtmst = [];
end
rmean = mean(rall,2);
mmean = mean(mall,2);
