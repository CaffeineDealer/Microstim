function [TC] = PDFalign(datapath,opt,type,flag,MS)





inf = load([datapath sprintf('N%sMS.mat',type)]); 
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
end
dir = (0:45:315);

TC = [];
for i = 1:size(MST,2)
    [rMT cMT] = max(MT(i).df);
    for j = 1:size(MST(i).df,2)
        [rMST cMST] = max(MST(i).df(:,j));
        idMT = cMST; %- cMT;
        for z = 1:size(MT(i).df,2)
%             TC(i).MST(j).bMT(:,z) = circshift(MT(i).df(:,z),idMT-1); % idMT(z) idMT-1
            dirnew = dir - dir(cMST);
            [rLB cLB] = find(dirnew < -180);
            ctlLB = size(dirnew(cLB),2);
            [rUB cUB] = find(dirnew > 135);
            ctlUB = size(dirnew(cUB),2);
            if ctlLB > 0
                TC(i).MST(j).bMT(:,z) = circshift(MT(i).df(:,z),-ctlLB);
%                 dirnew(cLB) = 45:45:ctlLB*(45);
%                 dirnew = circshift(dirnew,-ctlLB);
            elseif ctlUB > 0
                TC(i).MST(j).bMT(:,z) = circshift(MT(i).df(:,z),ctlUB);
%                 dirnew(cUB) = 45:45:ctlUB*(45);
%                 dirnew = circshift(dirnew,ctlUB);
            elseif ctlLB == 0 && ctlUB == 0
                TC(i).MST(j).bMT(:,z) = circshift(MT(i).df(:,z),0);
            end
        end
        TC(i).MST(j).MST = MST(i).df(:,j);
        TC(i).MST(j).MTmean = mean(TC(i).MST(j).bMT,2);
        TC(i).MST(j).MSTidx = cMST;     
%         TC(i).MST(j).dir = -((cMST - 1) * 45):45:((8 - cMST) * 45);
        TC(i).MST(j).dir = -180:45:135;
        TC(i).MST(j).err = std(TC(i).MST(j).bMT,[],2) / sqrt(size(TC(i).MST(j).bMT,2));
        TC(i).MST(j).fr = MST(i).FR(j,1);
    end
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

