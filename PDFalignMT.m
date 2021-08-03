function [TC] = PDFalignMT(datapath,opt,type,flag,MS)





inf = load([datapath sprintf('N%s.mat',type)]); 
MST = inf.info; 
clear info inf
inf = load([datapath sprintf('%s%s.mat',opt,type)]); 
MT = inf.info; 
clear info inf

if opt == 'S'
    inf = load([datapath sprintf('%s%s.mat','N',type)]);
    N = inf.info;
    for i = 1:size(N,2)
        for j = 1:size(N(i).xr,2)
            MT(i).df(:,j) = MT(i).xr(j).firing(:,N(i).pt(j));
        end
    end
end

if MS == 1
    MT(7) = [];
    MT(4) = []; 
    MT(1) = [];
    MST(7) = [];
    MST(4) = [];
    MST(1) = [];
end
dir = (0:45:315);
TC = [];
for i = 1:size(MST,2)
    for j = 1:size(MST(i).df,2)
        [rMST cMST] = max(MST(i).df(:,j));
        dirnew = dir - dir(cMST);
        [rLB cLB] = find(dirnew < -180);
        ctlLB = size(dirnew(cLB),2);
        [rUB cUB] = find(dirnew > 135);
        ctlUB = size(dirnew(cUB),2);
        if ctlLB > 0
            TC(i).MST.bMT(:,j) = circshift(MT(i).df(:,j),-ctlLB);
        elseif ctlUB > 0
            TC(i).MST.bMT(:,j) = circshift(MT(i).df(:,j),ctlUB);
        elseif ctlLB == 0 && ctlUB == 0
            TC(i).MST.bMT(:,j) = circshift(MT(i).df(:,j),0);
        end
        TC(i).MST.MSTidx(:,j) = cMST;
        TC(i).MST.fr(:,j) = MST(i).FR(j,1);
    end
%     TC(i).MST.MST = MST(i).df;
    TC(i).MST.MTmean = mean(TC(i).MST.bMT,2);
    TC(i).MST.dir = -180:45:135;
    TC(i).MST.err = std(TC(i).MST.bMT,[],2) / sqrt(size(TC(i).MST.bMT,2));
    if flag == 1
        if ceil(j/2) == 1
            r = 1; c = 2;
        else
            r = ceil(j/2);
            c = r;
        end
        figure
        for w = 1:j
            subplot(r,c,w)
            plot(TC(i).MST(w).dir,TC(i).MST(w).MST,'k','LineWidth',0.5)
            hold on
            plot(TC(i).MST(w).dir,TC(i).MST(w).MTmean,'b','LineWidth',1.5)
            xlim([min(TC(i).MST(w).dir) max(TC(i).MST(w).dir)])
        end
        title(sprintf('Mean MT TC^{%s}',type))
        xlabel 'Diff from MST PD^o'
        ylabel 'Firing Rate (Hz)'
    end
end