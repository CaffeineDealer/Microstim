function [x] = PosAtPDMST(datapath,type,MS,site)





% Function PosAtPD estimates the euclidean distance of a position from a
% reference point on the screen
%%%% Inputs %%%%
% x:    x-axis coordinates / eg. x = [-17;0;17;-17;0;17;-17;0;17];
% y:    y-axis coordinates / eg. y = [17;17;17;0;0;0;-17;-17;-17];
% ref:  reference coordinates / eg. ref = [17 17];
%%%% Outputs %%%%
% dist:   distance from reference point

% Written by Yavar Korkian on July.19.2021

switch type
    case 'trans'
        mot = 1;
    case 'spiral'
        mot = 2;
end
bin = 100;
Fs = 10000;

inf = load([datapath sprintf('N%sMS.mat',type)]);
MST = inf.info;
clear inf
inf = load([datapath sprintf('N%s.mat',type)]);
MTnms = inf.info;
clear inf
inf = load([datapath sprintf('S%s.mat',type)]);
MTms = inf.info;
clear inf

sppath = 'E:\MT_MST\SuperTuneSpkTrains\ms_rmvd\';
for i = 1:size(MTnms,2)
    for j = 1:size(MTnms(i).xr,2)
        MTnms(i).spike(j) = load([sppath MTnms(i).fname sprintf('N%d1spktrain.mat',MTnms(i).gch(j))]);
        MTms(i).spike(j) = load([sppath MTnms(i).fname sprintf('S%d1spktrain.mat',MTnms(i).gch(j))]);
    end
end

for i = 1:size(MTnms,2)
    for j = 1:size(MTnms(i).xr,2)
        MTms(i).df(:,j) = MTms(i).xr(j).firing(:,MTnms(i).pt(j));
    end
end

bl(:,1) = vertcat(MTnms.spnt);
bl(:,2) = vertcat(MTms.spnt);
bl = mean(bl(:));

if MS == 1
    MTnms(7) = [];
    MTnms(4) = [];
    MTnms(1) = [];
    MTms(7) = [];
    MTms(4) = [];
    MTms(1) = [];
end

NMS = nan(400,6);
MS = NMS;
z0 = 1; z1 = 1; z2 = 1; z3 = 1; z4 = 1; z5 = 1;
x = [];

for i = 1:size(MST,2)
    for z = 1:size(MST(i).df,2)
        FrMST = MST(i).xr(z).firing;
        [rval I] = max(FrMST(:));
        [pdr pdc] = ind2sub(size(FrMST),I);
        for j = 1:size(MTnms(i).xr,2)
            FrMTnms = MTnms(i).xr(j).firing;
            FrMTms = MTms(i).xr(j).firing;
            x(i).MST(z).dist(:,j) = pos2dist(pdc);
            x(i).MST(z).FrmMTnms(:,j) = FrMTnms(pdr,:)';
            x(i).MST(z).FrmMTms(:,j) = FrMTms(pdr,:)';
            % spike train
            SpMTnms = squeeze(MTnms(i).spike(j).spktrain(:,pdr,mot,:,:));
            SpMTms = squeeze(MTms(i).spike(j).spktrain(:,pdr,mot,:,:));
            % psth
            spbinNMS = []; spbinMS = [];
            for zz = 1:9 %
                spbinNMS(:,zz) = psth_sp(squeeze(SpMTnms(:,zz,:)),bin,Fs);
                spbinMS(:,zz) = psth_sp(squeeze(SpMTms(:,zz,:)),bin,Fs);
            end
            
            nms = x(i).MST(z).FrmMTnms(:,j);
            ms = x(i).MST(z).FrmMTms(:,j);
            b = x(i).MST(z).dist(:,j);
            upos = unique(b);
            for w = 1:size(upos,1)
                x(i).MST(z).pos(j).MT(w,1) = upos(w);
                x(i).MST(z).pos(j).MT(w,2) = mean(nms(b == upos(w))); % NMS
                x(i).MST(z).pos(j).MT(w,3) = mean(ms(b == upos(w))); % MS
                x(i).MST(z).pos(j).spNMS(:,w) = mean(spbinNMS(:,b == upos(w)),2);
                x(i).MST(z).pos(j).spMS(:,w) = mean(spbinMS(:,b == upos(w)),2);
                switch upos(w)
                    case 0
                        NMS(z0,1) = mean(nms(b == upos(w))); % NMS
                        MS(z0,1) = mean(ms(b == upos(w))); % MS
                        z0 = z0 + 1;
                    case 1
                        NMS(z1,2) = mean(nms(b == upos(w))); % NMS
                        MS(z1,2) = mean(ms(b == upos(w))); % MS
                        z1 = z1 + 1;
                    case 2
                        NMS(z2,3) = mean(nms(b == upos(w))); % NMS
                        MS(z2,3) = mean(ms(b == upos(w))); % MS
                        z2 = z2 + 1;
                    case 3
                        NMS(z3,4) = mean(nms(b == upos(w))); % NMS
                        MS(z3,4) = mean(ms(b == upos(w))); % MS
                        z3 = z3 + 1;
                    case 4
                        NMS(z4,5) = mean(nms(b == upos(w))); % NMS
                        MS(z4,5) = mean(ms(b == upos(w))); % MS
                        z4 = z4 + 1;
                    case 5
                        NMS(z5,6) = mean(nms(b == upos(w))); % NMS
                        MS(z5,6) = mean(ms(b == upos(w))); % MS
                        z5 = z5 + 1;
                end
            end
        end
    end
end

for w = 1:6
    a = []; b = [];
    c = 1;
    for i = 1:size(x,2)
        for z = 1:size(x(i).MST,2)
            for j = 1:size(x(i).MST(z).pos,2)
                if size(x(i).MST(z).pos(j).spNMS,2) >= w
                    a(:,c) = x(i).MST(z).pos(j).spNMS(1:24,w);
                    b(:,c) = x(i).MST(z).pos(j).spMS(1:24,w);
                    c = c + 1;
                end
            end
        end
    end
    sptnms(:,w) = mean(a,2);
    sptms(:,w) = mean(b,2);
end

for i = 1:size(NMS,2)
    n = NMS(~isnan(NMS(:,i)),i);
    s = MS(~isnan(MS(:,i)),i);
    errN(1,i) = (std(n)) / (sqrt(size(n,1))); 
    errS(1,i) = (std(s)) / (sqrt(size(s,1))); 
    n = []; s = [];
end

NMS = nanmean(NMS);
MS = nanmean(MS);

figure
dir = 1:6;
switch type
    case 'trans'
        col(1) = 'k';
        col(2) = 'g';
        tit = 'Translation';
    case 'spiral'
        col(1) = 'b';
        col(2) = 'r';
        tit = 'Spiral';
end
plot(dir,NMS,'Color',col(1),'LineWidth',2); hold on
plot(dir,MS,'Color',col(2),'LineWidth',2)

title(tit)
xlabel(sprintf('Position Deviation from %s',site))
ylabel 'Normalized Firing Rate'
legend({'NMS','MS'})
xticks([1 2 3 4 5 6])
xticklabels({'0','1','2','3','4','5'})

shade(dir,NMS+errN,sprintf('--%s',col(1)),dir,NMS-errN,sprintf('--%s',col(1)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(1))
shade(dir,MS+errS,sprintf(':%s',col(2)),dir,MS-errS,sprintf(':%s',col(2)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(2))

line([min(dir) max(dir)],[bl bl],'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--')

fig = figure;
time = 0:0.01:0.230;
c = 1;
for i = 1:6
    if i == 1
        subplot(6,6,1)
    else
        c = c + 7;
        subplot(6,6,c)
    end
    bar(time,sptnms(:,i),'BarWidth',1,'FaceColor',col(1),'EdgeColor',col(1)); hold on
    bar(time,sptms(:,i),'BarWidth',1,'FaceColor',col(2),'EdgeColor',col(2),'FaceAlpha',0.4)
    xlim([min(time) max(time)])
    gmin = min([min(sptnms(:,i)) min(sptms(:,i))]);
    gmax = max([max(sptnms(:,i)) max(sptms(:,i))]);
    ylim([0 gmax+1])
end
legend({'NMS','MS'})
han = axes(fig,'visible','off'); 
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
title(han,tit)
ylabel(han,'Firing Rate (spike/sec)')
xlabel(han,'Time after MS (sec)')
