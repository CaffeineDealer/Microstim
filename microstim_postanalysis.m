%% Start: Load N & S datasets 
clear all; clc
% fldname = 'E:\MT_MST\Microstim\Norm\MUA_SUA\';
fldname = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
param.mua_sua = 0;
param.type = 'trans';
param.option = '';
[N,S] = LoadDataset(fldname,param);
%% Vector average: Preferred Direction MT
svpath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\Norm\';
svpath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
[badCell,theta] = mtpd(N,S,svpath,1,param.option,param.type);
%% Vector average: Preferred Direction MST
svpath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\Norm\';
[badCell,theta] = mstpd(N,svpath,0,param.option,param.type);
%% RF MT
[Vx Vy diffV] = RF(S,N,badCell,type,theta,1);
% save(['E:\MT_MST\Microstim\PSTH\MUA_SUA\',sprintf('cmass_%s_%sn.mat',type(1),option)],'Vx','Vy','diffV')
%% RF MST
S = N;
[Vx Vy diffV] = RF(S,N,badCell,type,theta,0);
save(['E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\',sprintf('cmass_%s_%s.mat',type(1),option)],'Vx','Vy','diffV')
% save(['E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\',sprintf('mst_RFcmass_%s%s.mat',type,option)],'Vx','Vy','diffV')
%% Position @ PDirection
clear all; clc
datapath = 'E:\MT_MST\Microstim\Norm\MUA_SUA\';
type = 'spiral';
MS = 1;
MT = PosAtPD(datapath,type,MS,'MT');
MST = PosAtPDMST(datapath,type,MS,'MST');
%% MST Stimulus Response Function
clear all; clc
datapath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\'; %ms_rmvd_150300
type = 'trans';
MS = 1;
[rn mn] = StimRespMST(datapath,'N',type,MS);
[rs ms] = StimRespMST(datapath,'S',type,MS);

figure
plot(rn,smooth(rn,mn,0.2),'k','LineWidth',1.5); hold on
plot(rn,smooth(rn,ms,0.2),'g','LineWidth',1.5)
xlabel 'Worst --> Best'
ylabel 'Normalized MT Response'
title (type)
xlim([0 1])
legend({'NMS','MS'},'Location','North')
%%
clear all; clc
datapath = 'E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\';
MS = 0;
type = 'trans';
TCN = PDFalign(datapath,'N',type,0,MS);
TCS = PDFalign(datapath,'S',type,0,MS);
%% MT Tuning curves aligned with MST: Part I
clear all; clc
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\';
datapath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\';
MS = 1;
type = 'trans';
[Frmt(:,1) Frmst  err(:,1)] = FR_MTMST(datapath,'N',type,0,MS);
[Frmt(:,2),~,err(:,2)] = FR_MTMST(datapath,'S',type,0,MS);

Fr = [];
c = 1;
for i = 1:max(Frmst(:,2))
    z = Frmst(Frmst(:,2)==i,1);
    for j = 1:size(z,1)
        Fr(c,:) = [Frmt(i,:) z(j) err(i,:)];
        c = c + 1;
    end   
end
Fr = sortrows(Fr,3);

figure
errorbar(Fr(:,3),Fr(:,1),Fr(:,4),'k','LineWidth',1.5); hold on
errorbar(Fr(:,3),Fr(:,2),Fr(:,5),'b','LineWidth',1.5)
title(type)
legend({'NMS','MS'})
xlabel 'Normalized MST Firing Rate'
ylabel 'Normalized MT Firing Rate'
xlim([0 1])
ylim([0 1])
%% MT Tuning curves aligned with MST: Part II.I
dir = (-180:45:135);
for i = 1:size(TCN,2)
    if ceil(size(TCN(i).MST,2)/2) == 1
        r = 1; c = 2;
    else
        r = ceil(size(TCN(i).MST,2)/2);
        c = r;
    end
    figure
    for w = 1:size(TCN(i).MST,2)
        subplot(r,c,w)
        plot(dir,TCN(i).MST(w).MTmean,'b','LineWidth',1.5); %TCN(i).MST(w).dir
        hold on
        plot(dir,TCS(i).MST(w).MTmean,'r','LineWidth',1.5); % TCS(i).MST(w).dir
%         xlim([min(TCN(i).MST(w).dir) max(TCN(i).MST(w).dir)])
        xlim([-180 135])
        xlabel(sprintf('Deviation from MST^{#%d} PD^o',w))
    end
    title(sprintf('Mean MT TC^{%s}',type))
    ylabel 'Firing Rate (Hz)'
    legend({'NMS','MS'})
end
%% Normalize all MST tuning curves
clear all; clc
datapath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\';%ms_rmvd_150300
svpath = 'E:\MT_MST\Microstim\Norm\MUA_SUA\MST\';

MSTtnms = load([datapath sprintf('N%sMS.mat','trans')]); 
MSTsnms = load([datapath sprintf('N%sMS.mat','spiral')]);

normtype = 'all'; 

switch normtype
% Normalized by the overall maximum (motion type & NMS & MS all together)
    case 'all'
        for i = 1:size(MSTtnms.info,2)
            maxTS = [];
            maxTS(1,:) = max(MSTtnms.info(i).df);
            maxTS(2,:) = max(MSTsnms.info(i).df);
            x(i).gmax = max(maxTS)';
            for j = 1:size(MSTtnms.info(i).xr,2)
                MSTtnms.info(i).xr(j).firing = MSTtnms.info(i).xr(j).firing / x(i).gmax(j);
                MSTsnms.info(i).xr(j).firing = MSTsnms.info(i).xr(j).firing / x(i).gmax(j);
            end
            MSTtnms.info(i).df = MSTtnms.info(i).df ./ x(i).gmax';
            MSTsnms.info(i).df = MSTsnms.info(i).df ./ x(i).gmax';
        end
% Normalized by the overall maximum for each motion type (NMS & MS all together)
    case 'sep'
        for i = 1:size(MSTtnms.info,2)
            maxTS = [];
            maxTS(1,:) = max(MSTtnms.info(i).df);
            maxTS(2,:) = max(MSTsnms.info(i).df);
            for j = 1:size(MSTtnms.info(i).xr,2)
                MSTtnms.info(i).xr(j).firing = MSTtnms.info(i).xr(j).firing / maxTS(1,j);
                MSTsnms.info(i).xr(j).firing = MSTsnms.info(i).xr(j).firing / maxTS(2,j);
            end
            MSTtnms.info(i).df = MSTtnms.info(i).df ./ maxTS(1,:);
            MSTsnms.info(i).df = MSTsnms.info(i).df ./ maxTS(2,:);
        end
end

clear info
info = MSTtnms.info;
save([svpath sprintf('N%sMS.mat','trans')],'info')
clear info

info = MSTsnms.info;
save([svpath sprintf('N%sMS.mat','spiral')],'info')
clear info
%% MST Firing Rate: Translation vs. Spirals
clear all; clc
datapath = 'E:\MT_MST\Microstim\Norm\MUA_SUA\MST\';
inf = load([datapath sprintf('N%sMS.mat','trans')]);
MSTt = inf.info;
clear inf
inf = load([datapath sprintf('N%sMS.mat','spiral')]);
MSTs = inf.info;
clear inf

mstT = []; mstS = [];
for i = 1:size(MSTt,2)
    mstT = [mstT mean(MSTt(i).df)];
    mstS = [mstS mean(MSTs(i).df)];
end
offset = 0.05;
figure
plot(mstT,mstS,'k.','MarkerSize',20)
patch([0 1 0],[0 1 1],'b','FaceAlpha',.2)
patch([0 1 1],[0 1 0],'r','FaceAlpha',.2)
xlim([min([mstT mstS])-offset max([mstT mstS])+offset])
ylim([min([mstT mstS])-offset max([mstT mstS])+offset])
title 'Mean Firing Rate of MST cells'
xlabel 'Norm. Fr(Translation)'
ylabel 'Norm. Fr(Spirals)'
%% Normalize all MT tuning curves: Each cell across all stimulus
clear all; clc
datapath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\';%ms_rmvd_150300
svpath = 'E:\MT_MST\Microstim\Norm\MUA_SUA\';
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\';
% svpath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';

% load([datapath sprintf('N%sMS.mat',type)]); 
% save([svpath sprintf('N%sMS.mat',type)],'info')
% clear info

MTtnms = load([datapath sprintf('N%s.mat','trans')]);
MTtms = load([datapath sprintf('S%s.mat','trans')]);
MTsnms = load([datapath sprintf('N%s.mat','spiral')]);
MTsms = load([datapath sprintf('S%s.mat','spiral')]);

normtype = 'sep'; 

switch normtype
% Normalized by the overall maximum (motion type & NMS & MS all together)
    case 'all'
        for i = 1:size(MTtnms.info,2)
            maxTS = [];
            maxTS(1,:) = max(MTtnms.info(i).df);
            maxTS(2,:) = max(MTsnms.info(i).df);
            x(i).gmax = max(maxTS)';
            for j = 1:size(MTtnms.info(i).xr,2)
                MTtnms.info(i).xr(j).firing = MTtnms.info(i).xr(j).firing / x(i).gmax(j);
                MTtms.info(i).xr(j).firing = MTtms.info(i).xr(j).firing / x(i).gmax(j);
                MTsnms.info(i).xr(j).firing = MTsnms.info(i).xr(j).firing / x(i).gmax(j);
                MTsms.info(i).xr(j).firing = MTsms.info(i).xr(j).firing / x(i).gmax(j);
            end
            MTtnms.info(i).df = MTtnms.info(i).df ./ x(i).gmax';
            MTtms.info(i).df = MTtms.info(i).df ./ x(i).gmax';
            MTsnms.info(i).df = MTsnms.info(i).df ./ x(i).gmax';
            MTsms.info(i).df = MTsms.info(i).df ./ x(i).gmax';
        end
% Normalized by the overall maximum for each motion type (NMS & MS all together)
    case 'sep'
        for i = 1:size(MTtnms.info,2)
            maxTS = [];
            maxTS(1,:) = max(MTtnms.info(i).df);
            maxTS(2,:) = max(MTsnms.info(i).df);
            for j = 1:size(MTtnms.info(i).xr,2)
                MTtnms.info(i).xr(j).firing = MTtnms.info(i).xr(j).firing / maxTS(1,j);
                MTtms.info(i).xr(j).firing = MTtms.info(i).xr(j).firing / maxTS(1,j);
                MTsnms.info(i).xr(j).firing = MTsnms.info(i).xr(j).firing / maxTS(2,j);
                MTsms.info(i).xr(j).firing = MTsms.info(i).xr(j).firing / maxTS(2,j);
            end
            MTtnms.info(i).df = MTtnms.info(i).df ./ maxTS(1,:);
            MTtms.info(i).df = MTtms.info(i).df ./ maxTS(1,:);
            MTsnms.info(i).df = MTsnms.info(i).df ./ maxTS(2,:);
            MTsms.info(i).df = MTsms.info(i).df ./ maxTS(2,:);
            MTtnms.info(i).spnt = MTtnms.info(i).ctl(:,2) ./ maxTS(1,:)';
            MTtms.info(i).spnt = MTtms.info(i).ctl(:,2) ./ maxTS(1,:)';
            MTsnms.info(i).spnt = MTsnms.info(i).ctl(:,2) ./ maxTS(1,:)';
            MTsms.info(i).spnt = MTsms.info(i).ctl(:,2) ./ maxTS(1,:)';
        end
end

clear info
info = MTtnms.info;
save([svpath sprintf('N%s.mat','trans')],'info')
clear info
info = MTtms.info;
save([svpath sprintf('S%s.mat','trans')],'info')
clear info
info = MTsnms.info;
save([svpath sprintf('N%s.mat','spiral')],'info')
clear info
info = MTsms.info;
save([svpath sprintf('S%s.mat','spiral')],'info')
clear info
%% PD difference between MST(NoMS) & MT(NoMS & MS): Part II
clear all; clc
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\Norm\'; % MUA-MST-V1V2 / ms_rmvd_150300 / MUA_SUA
datapath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\';
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA_V1\';
[pdDiff Fr] = prefDiff(datapath,20);
%% MT Tuning curves aligned with MST: Part II.II
clear all; clc
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\';%ms_rmvd_150300
datapath = 'E:\MT_MST\Microstim\Norm\MUA_SUA\';
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
site = 'MST';
type = 'trans';
MS = 1;

switch site
    case 'MST'
        TCN = PDFalign(datapath,'N',type,0,MS);
        TCS = PDFalign(datapath,'S',type,0,MS);
    case 'MT'
        TCN = PDFalignMT(datapath,'N',type,0,MS);
        TCS = PDFalignMT(datapath,'S',type,0,MS);
end

[N S errN errS dir] = MTalignedTC(TCN,TCS,type,site);
%%
for i = 1:size(N,2)
    for j = 1:size(N(i).xr,2)
        S(i).df(:,j) = S(i).xr(j).firing(:,N(i).pt(j));
    end
end

%% Master Plot(Fr & PD): MT vs MST 
clear all; clc
datapath = 'E:\MT_MST\Microstim\Norm\MUA_SUA\';
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
type = 'trans';
ty = 'spiral';
MS = 1;
inf = load([datapath sprintf('N%sMS.mat',type)]); 
MSTt = inf.info; 
clear info inf
inf = load([datapath sprintf('N%s.mat',type)]); 
MTtnms = inf.info; 
clear info inf
inf = load([datapath sprintf('S%s.mat',type)]); 
MTtms = inf.info; 
clear info inf
%
inf = load([datapath sprintf('N%sMS.mat',ty)]); 
MSTs = inf.info; 
clear info inf
inf = load([datapath sprintf('N%s.mat',ty)]); 
MTsnms = inf.info; 
clear info inf
inf = load([datapath sprintf('S%s.mat',ty)]); 
MTsms = inf.info; 
clear info inf

for i = 1:size(MTtnms,2)
    for j = 1:size(MTtnms(i).xr,2)
        MTtms(i).df(:,j) = MTtms(i).xr(j).firing(:,MTtnms(i).pt(j));
        MTsms(i).df(:,j) = MTsms(i).xr(j).firing(:,MTsnms(i).pt(j));
    end
end

if MS == 1
    MTtnms(7) = [];
    MTtnms(4) = [];
    MTtnms(1) = [];
    MTtms(7) = [];
    MTtms(4) = [];
    MTtms(1) = [];
    MTsnms(7) = [];
    MTsnms(4) = [];
    MTsnms(1) = [];
    MTsms(7) = [];
    MTsms(4) = [];
    MTsms(1) = [];
end

[cort] = masterPlot(MSTt,MTtnms,MTtms,'Translation'); %x
[cors] = masterPlot(MSTs,MTsnms,MTsms,'Spiral');
%% 
clear all; clc
datapath = 'E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\';
load([datapath,'spiralMST.mat'])
load([datapath,'spiralMT.mat'])
Nmst = mean(Nmst,2);
Nmt = mean(Nmt,2);
Smt = mean(Smt,2);
Smst = mean(Smst,2);

N = Nmst - Nmt;
S = Smst - Smt;
NMS = diag(N)*100;
MS = diag(S)*100;
type = 'Spirals';
gmax = ceil(max(max(NMS(:)),max(MS(:))))+180;
gmin = floor(min(min(NMS(:)),min(MS(:))));

figure
x = (0:45:180);
for i = 1:5
    plot(x(i)+NMS(i,i),x(i)+NMS(i,i),'b.','MarkerSize',40); hold on
    plot(x(i)+MS(i,i),x(i)+MS(i,i),'r.','MarkerSize',40); hold on
end
% refline(1,0)
grid on
xlim([gmin-5 gmax+5])
ylim([gmin-5 gmax+5])
title(sprintf('%s',type))
legend({'NMS','MS'},'Location','NorthWest')
xlabel 'Deviation from MT PD^{o}'
ylabel 'Deviation from MST PD^{o}'
set(gca,'XTick',[0 45 90 135 180],'YTick',[0 45 90 135 180])
%% Compute & Plot PSTH
sst = horzcat(S.sp);
nst = horzcat(N.sp);
bin = 100;
spike = [];
for j = 1:size(sst,2)
    Ssp = []; Ssp_bl = []; Nsp = []; Nsp_bl = [];
    for i = 1:4
        Ssp(i,:) = psth_sp(sst(j).spktrain,bin,i);
        Ssp_bl(i,:) = psth_sp(sst(j).spktrain_bl,bin,i);
        Nsp(i,:) = psth_sp(nst(j).spktrain,bin,i);
        Nsp_bl(i,:) = psth_sp(nst(j).spktrain_bl,bin,i);
    end
    spike(j).Ssp = mean(Ssp);
    spike(j).Ssp_bl = mean(Ssp_bl);
    spike(j).Nsp = mean(Nsp);
    spike(j).Nsp_bl = mean(Nsp_bl);
    sizeStim(j,1) = size(spike(j).Ssp,2);
    sizeBl(j,1) = size(spike(j).Ssp_bl,2);
end
ctl(1,1) = min(sizeBl);
ctl(1,2) = min(sizeStim);
for j = 1:size(sst,2)
    spike(j).Ssp_bl = spike(j).Ssp_bl(1,sizeBl(j,1)-ctl(1,1)+1:end);
    spike(j).Nsp_bl = spike(j).Nsp_bl(1,sizeBl(j,1)-ctl(1,1)+1:end);
    spike(j).Ssp = spike(j).Ssp(1,1:ctl(1,2));
    spike(j).Nsp = spike(j).Nsp(1,1:ctl(1,2)); 
end
Ssp = vertcat(spike.Ssp);
Ssp_bl = vertcat(spike.Ssp_bl);
a = zeros(size(Ssp,1),14);
S = [Ssp_bl a Ssp];
Nsp = vertcat(spike.Nsp);
Nsp_bl = vertcat(spike.Nsp_bl);
N = [Nsp_bl a Nsp];
Smean = mean(S);
Nmean = mean(N);
t = (-ctl(1)+1)*10:10:(ctl(2)+14)*10;
figure
plot(t,smooth(Nmean,3),'LineWidth',1,'Color','k')
hold on
plot(t,smooth(Smean,3),'LineWidth',1,'Color','g')
title(sprintf('MT-%s',type))
xlim([-400 400])
xlabel('Time(ms)')
ylabel('Firing Rate')
legend({'No-Stim' 'Stim'})
%% Max firing rate: No-MicroStim vs Microstim
npd = horzcat(N.df);
spdgmax = horzcat(S.df);
% pos = vertcat(N.pt);
% dir = vertcat(N.Dir); dir(:,2) = [];
n = 1;
for i = 1:size(N,2)
    for j = 1:size(S(i).xr,2)
        spd(:,n) = S(i).xr(j).firing(:,N(i).pt(j));
        n = n + 1;
    end
end
ctl = [];
for i = 1:size(spd,2) 
    if spd(:,i) == spdgmax(:,i) 
        ctl(1,i) = 1; 
    else
        ctl(1,i) = 0; 
    end
end
% cmaxmap = zeros(size(spd,2),3);
% cmaxmap(ctl==0,1) = 1; 
% scatter(max(spd(:,ctl==0)),max(spdgmax(:,ctl==0)),25,'r','filled'); refline(1,0)

[frm frpos] = max(npd);
for i = 1:size(frpos,2)
    frm(2,i) = spd(frpos(i),i); % Firing rate for the same position
end
frm(3,:) = max(spd); % Max within 8 dir
frm(4,:) = 0;
frm(4,frm(2,:)==frm(3,:)) = 1;
cmap = zeros(size(frm,2),3);
figure
switch type
    case 'trans'
        cmap(:,1) = 0;
        cmap(:,2) = 0;
        cmap(:,3) = 1;
        cmap(frm(4,:)==1,2) = 1;
        cmap(frm(4,:)==1,3) = 0;
        subplot(1,3,2)
        scatter(frm(1,:),frm(3,:),25,'b','filled'); refline(1,0)
        title 'Local Maxima'
%         xlim([0 210]);
%         ylim([0 210]);
        xlabel 'Fr^{NMS} (spk/sec)'
        %         scatter(frm(1,:),frm(2,:),'MarkerEdgeColor','g','MarkerFaceColor','b'); refline(1,0)
    case 'spiral'
        cmap(frm(4,:)==1,1) = 1;
        subplot(1,3,2)
        scatter(frm(1,:),frm(3,:),25,'k','filled'); refline(1,0)
        title 'Local Maxima'
%         xlim([0 210]);
%         ylim([0 210]);
        xlabel 'Fr^{NMS} (spk/sec)'
%         scatter(frm(1,:),frm(2,:),'MarkerEdgeColor','r','MarkerFaceColor','k'); refline(1,0)
end
subplot(1,3,1)
scatter(frm(1,:),frm(2,:),25,cmap,'filled'); refline(1,0)
title(sprintf('Resp^{%s}_{MT} P2P',type))
% xlim([0 210])
% ylim([0 210])
ylabel 'Fr^{MS} (spk/sec)'
subplot(1,3,3)
scatter(max(spd(:,ctl==0)),max(spdgmax(:,ctl==0)),25,'m','filled'); refline(1,0)
title 'Global Maxima'
% title(sprintf('Resp^{%s}_{MT}',type))
% xlim([0 210]); xlabel 'Fr^{NMS} (spk/sec)' % max(frm(1,:))
% ylim([0 210]); ylabel 'Fr^{MS} (spk/sec)' % max(frm(2,:))
%% PD difference between MST(NoMS) & MT(NoMS & MS): Part I
switch option
    case 'MS'
        tr = struct;
        trMST = zeros(max(theta(:,7)),2);
        for i = 1:max(theta(:,7))
            a = find(theta(:,7)==i);
            tr(i).thetaFr(:,1:2) = theta(a,1:2); % PD for MS & NoMS
            tr(i).thetaFr(:,3:4) = theta(a,5:6); % Firing rate for MS & NoMS
            tr(i).firingS = spd(:,a);
            tr(i).firingN = npd(:,a);
            clear a
            [trMST(i,1),~,~] = vec_avg(tr(i).thetaFr(:,3)',tr(i).thetaFr(:,1)'); % MS
            [trMST(i,2),~,~] = vec_avg(tr(i).thetaFr(:,4)',tr(i).thetaFr(:,2)'); % NoMS
        end
        save(['E:\MT_MST\Microstim\PSTH\msr\mt-mst\',sprintf('trMST_%s_%s.mat',type(1),option)],'trMST','tr')        
end
%% PD difference between MST(NoMS) & MT(NoMS & MS): Part I - Updated
switch option
    case 'MS'
        tr = struct;
        trMST = zeros(max(theta(:,4)),1);
        for i = 1:max(theta(:,4))
            a = find(theta(:,4)==i);
            tr(i).thetaFr(:,1) = theta(a,1); % PD for NoMS
            tr(i).thetaFr(:,2) = theta(a,3); % Firing rate for NoMS
            tr(i).firingN = npd(:,a);
            clear a
            [trMST(i,1),~,~] = vec_avg(tr(i).thetaFr(:,2)',tr(i).thetaFr(:,1)'); % NoMS
        end
        save(['E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\',sprintf('trMST_%s_%s.mat',type(1),option)],'trMST','tr')      % MUA-MST-V1V2  
end
%% PD difference between MST(NoMS) & MT(NoMS & MS): Part II
clear all; clc
datapath = 'E:\MT_MST\Microstim\PSTH\MUA_SUA\Norm\'; % MUA-MST-V1V2 / ms_rmvd_150300 / MUA_SUA
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\';
% datapath = 'E:\MT_MST\Microstim\PSTH\MUA_V1\';
[pdDiff Fr] = prefDiff(datapath,20);
%% PD difference between MST(NoMS) & MT(NoMS & MS): Part II
clear all; clc
% load(['E:\MT_MST\Microstim\PSTH\MUA_SUA\' 'trMST_t_MS.mat']); trMSTt = trMST; clear trMST % MUA-MST-V1V2
% load(['E:\MT_MST\Microstim\PSTH\MUA_SUA\' 'trMST_s_MS.mat']); trMSTs = trMST; clear trMST
load(['E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\' 'theta_t_MS.mat']); trMSTt = theta; clear trMST % MUA-MST-V1V2
load(['E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\' 'theta_s_MS.mat']); trMSTs = theta; clear trMST
load(['E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\' 'theta_t_.mat']); trMTt = theta; clear theta
load(['E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\' 'theta_s_.mat']); trMTs = theta; clear theta

pdDiff = zeros(size(trMTt,1),4);
for i = 1:max(trMTt(:,7))
    mtidx = find(trMTt(:,7)==i);
    mstT = ones(size(mtidx,1),1) * trMSTt(i,1); % trMSTt(i,2)
    mstS = ones(size(mtidx,1),1) * trMSTs(i,1); % trMSTs(i,2)
    pdDiff(mtidx,1) = abs(angdiff(mstT,trMTt(mtidx,1))); % MST Trans NoMS vs. MT Trans MS
    pdDiff(mtidx,2) = abs(angdiff(mstT,trMTt(mtidx,2))); % MST Trans NoMS vs. MT Trans NoMS
    pdDiff(mtidx,3) = abs(angdiff(mstS,trMTs(mtidx,1))); % MST Spiral NoMS vs. MT Spiral MS
    pdDiff(mtidx,4) = abs(angdiff(mstS,trMTs(mtidx,2))); % MST Spiral NoMS vs. MT Spiral NoMS
end
pdDiff = rad2deg(pdDiff);
pdD(:,1) = rad2deg(trMTt(:,3));
pdD(:,2) = rad2deg(trMTs(:,3));
% save('E:\MT_MST\Microstim\PSTH\MUA_V1\pdDiff.mat','pdDiff')
% save('E:\MT_MST\Microstim\PSTH\MUA_V1\pdD.mat','pdD')
%%%% Plot Section %%%%
figure
subplot(1,2,1)
scatter(pdDiff(:,2),pdDiff(:,1),'MarkerEdgeColor','g','MarkerFaceColor','b'); refline(1,0)
xlim([0 181]); ylim([0 181])
title('PD_{Diff}^{Trans}')
xlabel('(MST_{NoMS} - MT_{NoMS})^o')
ylabel('(MST_{NoMS} - MT_{MS})^o')
subplot(1,2,2)
scatter(pdDiff(:,4),pdDiff(:,3),'MarkerEdgeColor','r','MarkerFaceColor','k'); refline(1,0)
xlim([0 181]); ylim([0 181])
title('PD_{Diff}^{Spiral}')
xlabel('(MST_{NoMS} - MT_{NoMS})^o')
%% Histogram of PD difference for translation and spiral
clear all; clc
cd E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\
option = '';
load(sprintf('theta_t_%s.mat',option)); theta_T = theta; clear theta
load(sprintf('theta_s_%s.mat',option)); theta_S = theta; clear theta
ctlxmin = min(min([theta_S(:,5) theta_T(:,5)]));
ctlxmax = max(max([theta_S(:,5) theta_T(:,5)]));
ctlymin = min(min([theta_S(:,6) theta_T(:,6)]));
ctlymax = max(max([theta_S(:,6) theta_T(:,6)]));

figure
subplot(2,2,1)
polarhistogram(theta_T(:,3),20,'FaceColor','b','EdgeColor','g','FaceAlpha',0.7,'EdgeAlpha',.5)
title 'P^{Trans}_{Dir}(S-N)';
subplot(2,2,2)
polarhistogram(theta_S(:,3),20,'FaceColor','k','EdgeColor','r','FaceAlpha',0.7,'EdgeAlpha',.5)
title 'P^{Spiral}_{Dir}(S-N)';
subplot(2,2,3)
scatter(theta_T(:,6),theta_T(:,5),'MarkerEdgeColor','g','MarkerFaceColor','b'); refline(1,0)
xlim([0 max(theta_T(:,5))+5])
ylim([0 max(theta_T(:,6))+5])
title 'R^{Trans}_{peak}'
ylabel 'MicroStim(spk/sec)' 
xlabel 'No MicroStim(spk/sec)'
subplot(2,2,4)
scatter(theta_S(:,6),theta_S(:,5),'MarkerEdgeColor','r','MarkerFaceColor','k'); refline(1,0)
% xlim([0 max(theta_S(:,5))+5])
% ylim([0 max(theta_S(:,6))+5])
xlim([0 max(theta_T(:,5))+5])
ylim([0 max(theta_T(:,6))+5])
title 'R^{Spiral}_{peak}'
ylabel 'MicroStim(spk/sec)' 
xlabel 'No MicroStim(spk/sec)'
%% Normalized Response as a ss of Position - based on max firing rate - RF Center of Mass - RF Shift
% Part I: Max firing rate btw 8 dir @ each 9 pos
% Part II: Center of Mass
% Part III: RF Center of Mass
% Part IV: RF Shift

% Max firing rate between 8 directions @ each 9 position
% S = N;
[Vx Vy diffV] = RF(S,N,badCell,type,theta,0);
save(['E:\MT_MST\Microstim\PSTH\MUA_V1\',sprintf('mst_RFcmass_%s%s.mat',type,option)],'Vx','Vy','diffV')
%%
[Vx Vy diffV] = RFcentered(S,N,badCell,type,theta,1);
%% Tuning curves: MT & MST : Translation & Spiral : No-MicroStim & MicroStim
clear all; clc
load('E:\MT_MST\Microstim\PSTH\MUA_SUA\ctl_MT.mat'); % msr\mt-mst
load('E:\MT_MST\Microstim\PSTH\MUA_SUA\MT.mat')
load('E:\MT_MST\Microstim\PSTH\MUA_SUA\pdDiff.mat')
load('E:\MT_MST\Microstim\PSTH\MUA_SUA\pdD.mat')
load('E:\MT_MST\Microstim\PSTH\MUA_SUA\trMST_t_MS.mat','tr'); trt = tr; clear tr
load('E:\MT_MST\Microstim\PSTH\MUA_SUA\trMST_s_MS.mat','tr'); trs = tr; clear tr
type = 'spiral';
names = {'ytu312b','ytu316a','ytu323a','ytu329a','ytu335a','ytu337a'};
dirs = deg2rad(0:45:315);
info = struct;
for i = 1:size(names,2) 
    info(i).MT = find(ctl_MT(:,2)==i);
    info(i).PDdiffT = pdDiff(info(i).MT,2); % MST(NoMS) - MT(NoMS) for Translation
    info(i).PDdiffS = pdDiff(info(i).MT,4); % MST(NoMS) - MT(NoMS) for Spiral
    info(i).PDdT = pdD(info(i).MT,1); % MT(NoMS) - MT(MS) for Translation
    info(i).PDdS = pdD(info(i).MT,2); % MT(NoMS) - MT(MS) for Spiral
    info(i).MTcell = load(['E:\MT_MST\Microstim\Cell_Lib\MST_MU_MT_SU\MT\','CLib_',names{1,i}(1:end-1),'.mat']); % MS_removed\Candidate Cells
    info(i).MSTcell = load(['E:\MT_MST\Microstim\Cell_Lib\MST_MU_MT_SU\MST\','CLib_',names{1,i}(1:end-1),'.mat']); % MS_removed\MS Cells
    count(i,1) = size(info(i).MT,1);
    count(i,2) = size(trt(i).firingN,2);
    count(i,3) = count(i,1) + count(i,2);
end
%%
% Plot Section
for i = 1:size(names,2)
    figure
    p = ceil(sqrt(count(i,3)));
    ct = 0;
    for j = 1:count(i,2)
        ct = ct + 1;
        subplot(p,p,ct)
        switch type
            case 'trans'
                plot(rad2deg(dirs),trt(i).firingN(:,j),'Color','m','LineWidth',1,'LineStyle','--'); hold on
            case 'spiral'
                plot(rad2deg(dirs),trs(i).firingN(:,j),'Color','m','LineWidth',1,'LineStyle','--'); hold on
        end
        xlim([0 315])
        title(sprintf('MST-%d',info(i).MSTcell.CLib(j)))
    end
    for z = 1:count(i,1)
        ct = ct + 1;
        subplot(p,p,ct)
        switch type
            case 'trans'
                plot(rad2deg(dirs),MT.ri_NMS_T(:,info(i).MT(z)),'Color','b','LineWidth',1,'LineStyle','--'); hold on
                plot(rad2deg(dirs),MT.ri_MS_T(:,info(i).MT(z)),'Color','b','LineWidth',2)
%                 title(sprintf('{MT_{N}- MT_{S}} = %.f^o / {MST_{N}- MT_{N}} = %.f^o',info(i).PDdT(z,1),info(i).PDdiffT(z,1)))
            case 'spiral'
                plot(rad2deg(dirs),MT.ri_NMS_S(:,info(i).MT(z)),'Color','k','LineWidth',1,'LineStyle','--'); hold on
                plot(rad2deg(dirs),MT.ri_MS_S(:,info(i).MT(z)),'Color','k','LineWidth',2)
%                 title(sprintf('{MT_{N}- MT_{S}} = %.f^o / {MST_{N}- MT_{N}} = %.f^o',info(i).PDdS(z,1),info(i).PDdiffS(z,1)))
        end
        xlim([0 315])
        title(sprintf('MT-%d',info(i).MTcell.CLib(z)))
    end
    xlabel 'Direction'; ylabel 'Fr (spk/sec)';
    suptitle(sprintf('Session %d',i))
end
%% Average MST cells
% Plot Section
count(:,2) = 1;
for i = 1:size(names,2)
    figure
    p = ceil(sqrt(count(i,3)));
    ct = 0;
    for j = 1:count(i,2)
        ct = ct + 1;
        subplot(p,p,ct)
        switch type
            case 'trans'
                plot(rad2deg(dirs),mean(trt(i).firingN,2),'Color','m','LineWidth',1,'LineStyle','--'); hold on
            case 'spiral'
                plot(rad2deg(dirs),mean(trs(i).firingN,2),'Color','m','LineWidth',1,'LineStyle','--'); hold on
        end
        xlim([0 315])
        title(sprintf('MST-%d',info(i).MSTcell.CLib(j)))
    end
    for z = 1:count(i,1)
        ct = ct + 1;
        subplot(p,p,ct)
        switch type
            case 'trans'
                plot(rad2deg(dirs),MT.ri_NMS_T(:,info(i).MT(z)),'Color','b','LineWidth',1,'LineStyle','--'); hold on
                plot(rad2deg(dirs),MT.ri_MS_T(:,info(i).MT(z)),'Color','b','LineWidth',2)
%                 title(sprintf('{MT_{N}- MT_{S}} = %.f^o / {MST_{N}- MT_{N}} = %.f^o',info(i).PDdT(z,1),info(i).PDdiffT(z,1)))
            case 'spiral'
                plot(rad2deg(dirs),MT.ri_NMS_S(:,info(i).MT(z)),'Color','k','LineWidth',1,'LineStyle','--'); hold on
                plot(rad2deg(dirs),MT.ri_MS_S(:,info(i).MT(z)),'Color','k','LineWidth',2)
%                 title(sprintf('{MT_{N}- MT_{S}} = %.f^o / {MST_{N}- MT_{N}} = %.f^o',info(i).PDdS(z,1),info(i).PDdiffS(z,1)))
        end
        xlim([0 315])
        title(sprintf('MT-%d',info(i).MTcell.CLib(z)))
    end
    xlabel 'Direction'; ylabel 'Fr (spk/sec)';
    suptitle(sprintf('Session %d',i))
end
%% Log for discarding bad cells
pathInfo = 'E:\MT_MST\Microstim\PSTH\n\';
dislog(1).trans = loadInfo(pathInfo,'trans',''); % MT trans
dislog(1).spiral = loadInfo(pathInfo,'spiral',''); % MT spiral
dislog(2).trans = loadInfo(pathInfo,'trans','MS'); % MST trans
dislog(2).spiral = loadInfo(pathInfo,'spiral','MS'); % MST spiral
for i = 1:7
    dislog(1).c(i).c = horzcat(dislog(1).trans(i).log,dislog(1).spiral(i).log(:,2:3)); % MT c1: Ch c2&c3: trans c4&c5: spiral
    dislog(2).c(i).c = horzcat(dislog(2).trans(i).log,dislog(2).spiral(i).log(:,2:3)); % MST c1: Ch c2&c3: trans c4&c5: spiral
end
%% VonMises Distribution 
dirs = linspace(-pi,pi,8+1);
for i = 1:size(theta_T,1)
    ri_MS_T(:,i) = vonMises([0.5 theta_T(i,1) theta_T(i,5) 0],dirs);
    ri_NMS_T(:,i) = vonMises([0.5 theta_T(i,2) theta_T(i,6) 0],dirs);
    ri_MS_S(:,i) = vonMises([0.5 theta_S(i,1) theta_S(i,5) 0],dirs);
    ri_NMS_S(:,i) = vonMises([0.5 theta_S(i,2) theta_S(i,6) 0],dirs);
end
%% Save ch & pos 
criteriaPD = 1; % 1: Change in PD   0: No change in PD
cond = 'S'; 
switch criteriaPD
    case 1
        flag = 'Y';
    case 0 
        flag = 'N';
end   
switch cond
    case 'S'
        criteriaPos = 4;
    case 'N'
        criteriaPos = 5;
end
for i = 1:size(names,2)
    chunum(:,1) = ctl(ctl(:,2)==i & ctl(:,3)==criteriaPD,1);
    chunum(:,2) = 1;
    pos = ctl(ctl(:,2)==i & ctl(:,3)==criteriaPD,criteriaPos);
    save(['E:\MT_MST\Microstim\PSTH\',sprintf('%sChunum%s%s.mat',names{i}(1:end-1)),cond,flag],'chunum')
    save(['E:\MT_MST\Microstim\PSTH\',sprintf('%sPos%s%s.mat',names{i}(1:end-1),cond,flag)],'pos')
    clear chunum pos
end
%% Create MT.mat
clear all; clc
cd E:\MT_MST\Microstim\PSTH\MUA_SUA
name = {'mt_trans.mat','mt_spiral.mat'};
dirs = deg2rad(0:45:315);
for z = 1:size(name,2)
    load(name{z})
    yMS = spd;
    yNMS = npd;
    for i = 1:size(theta,1)
        x(z).riMS(:,i) = yMS(:,i);
        x(z).riNMS(:,i) = yNMS(:,i);
    end
    x(z).fname = name{z};
    clear theta spd npd yMS yNMS
end
MT.ri_MS_T = x(1).riMS;
MT.ri_NMS_T = x(1).riNMS;
MT.ri_MS_S = x(2).riMS;
MT.ri_NMS_S = x(2).riNMS;
save(['E:\MT_MST\Microstim\PSTH\MUA_SUA\MT.mat'],'MT')
%% CMass stuff
VxMT = Vx_s_mt;
VyMT = Vy_s_mt;

VxMST = Vx_s_mst;
VyMST = Vy_s_mst;

c = 1;
for i = 1:max(ctl_mst)
    mtid = find(ctl_mt(:,1)==i);
    mstid = find(ctl_mst(:,1)==i);
    for j = 1:size(mstid,1)
        x2 = VxMST(mstid(j));
        y2 = VyMST(mstid(j));
        for z = 1:size(mtid,1)
            xS = VxMT(mtid(z),1);
            yS = VyMT(mtid(z),1);
            dt(c,1) = Edist([xS x2],[yS y2]);
            xN = VxMT(mtid(z),2);
            yN = VyMT(mtid(z),2);
            dt(c,2) = Edist([xN x2],[yN y2]);
            c = c + 1;
        end
    end
end
% dt_trans = dt;
dt_spiral = dt;
%%
[ht,pt] = kstest2(dt_trans(:,1),dt_trans(:,2));
[hs,ps] = kstest2(dt_spiral(:,1),dt_spiral(:,2));

figure
subplot(2,2,1)
F1 = cdfplot(dt_trans(:,1));
set(F1,'LineWidth',2,'Color','b')
hold on
F2 = cdfplot(dt_trans(:,2));
set(F2,'LineWidth',2,'Color','g')
legend({'MS','NMS'},'Location','East')
grid off
title 'CDF Center of Mass'
xlabel ''
ylabel 'Probability'

subplot(2,2,2)
F1 = cdfplot(dt_spiral(:,1));
set(F1,'LineWidth',2,'Color','k')
hold on
F2 = cdfplot(dt_spiral(:,2));
set(F2,'LineWidth',2,'Color','r')
legend({'MS','NMS'},'Location','East')
grid off
title 'CDF Center of Mass'
xlabel ''
ylabel ''

bin = 20;
subplot(2,2,3)
histogram(dt_trans(:,1),bin,'FaceColor','b')
hold on
histogram(dt_trans(:,2),bin,'FaceColor','g')
title(sprintf('Dev. from MST^{CM}_{Trans}  p-val = %.2f',pt))
xlabel('Distance')
ylabel('#of pairs')
% legend({'MS','NMS'})
xlim([0 max(dt_trans(:))+1])

subplot(2,2,4)
histogram(dt_spiral(:,1),bin,'FaceColor','k')
hold on
histogram(dt_spiral(:,2),bin,'FaceColor','r')
title(sprintf('Dev. from MST^{CM}_{Spiral}  p-val = %.2f',ps))
% legend({'MS','NMS'})
xlim([0 max(dt_spiral(:))+1])
%%
dt = vertcat(dt_trans,dt_spiral);
[h,p] = kstest2(dt(:,1),dt(:,2));
bin = 50;

figure
subplot(2,1,1)
F1 = cdfplot(dt(:,1));
set(F1,'LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
F2 = cdfplot(dt(:,2));
set(F2,'LineWidth',2,'Color',[0.4660 0.6740 0.1880])
legend({'MS','NMS'},'Location','East')
grid off
title 'CDF Center of Mass Combined'
xlabel ''
ylabel 'Probability'

subplot(2,1,2)
histogram(dt(:,1),bin,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.8)
hold on
histogram(dt(:,2),bin,'FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.8)
title(sprintf('Deviation from MST^{CM}_{Combined}  p-val = %.2f',p))
xlabel 'Distance'
ylabel '#of pairs'
% legend({'MS','NMS'})
xlim([0 max(dt(:))+1])
%%
dtdiff(:,1) = dt_trans(:,2) - dt_trans(:,1);
dtdiff(:,2) = dt_spiral(:,2) - dt_spiral(:,1);
[h,p] = kstest2(dtdiff(:,1),dtdiff(:,2));

bin = 50;
figure
subplot(2,1,1)
F1 = cdfplot(dtdiff(:,1));
set(F1,'LineWidth',2,'Color','b')
hold on
F2 = cdfplot(dtdiff(:,2));
set(F2,'LineWidth',2,'Color','k')
legend('Trans','Spiral','Location','East')
grid off
title 'CDF Center of Mass'
xlabel ''
ylabel 'Probability'

subplot(2,1,2)
histogram(dtdiff(:,1),bin,'FaceColor','b','FaceAlpha',0.6)
hold on
histogram(dtdiff(:,2),bin,'FaceColor','k','FaceAlpha',0.8)
% legend('Trans','Spiral')
title(sprintf('Distance of Diff from MST^{CM}  p-val = %.2f',p))
xlabel('Distance')
ylabel('#of pairs')

%%
dt = vertcat(dtdiff(:,1),dtdiff(:,2));
figure
histogram(dt,bin,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1)
title 'Distance of Diff from MST^{CM}_{Combined}'
xlabel 'Distance'
ylabel '#of pairs'
%% Fixation Experiment 
clear all; clc
names = {'ytu310a','ytu316a','ytu321a','ytu329a','ytu333a','ytu335a','ytu337a'};
frn = []; frs = [];
for i = 1:length(names)
    fn = names{i};
    a = struct2array(load(['E:\MT_MST\SuperTuneFiringMatrix\fixation\',fn(1:end-1),'cfr.mat']));
    frn = [frn;a(:,:,1)]; %N
    frs = [frs;a(:,:,2)]; %S
    ncell(i,1) = size(a,1);
end

figure
bl = mean([frn frs(:,2)],2);
ms = frs(:,1);
scatter(bl,ms,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor','k')
refline(1,0)
title(sprintf('MT Firing Rate NoStim (n=%d)',size(frs,1)))
xlabel 'Baseline (spk/sec)'
ylabel 'Microstim (spk/sec)'
xlim([0 25])
ylim([0 25])

%% Fixation vs. PD
datapath = 'E:\MT_MST\Microstim\Norm\MUA_SUA\';
type = 'trans';
inf = load([datapath sprintf('N%s.mat',type)]);
MTnms = inf.info;
clear inf
inf = load([datapath sprintf('S%s.mat',type)]);
MTms = inf.info;
clear inf

MTnms(5) = []; MTms(5) = [];
MTnms(2) = []; MTms(2) = [];

DirNMS = cell2mat({MTnms(:).Dir}');
Dir = indx2dir(DirNMS(:,1));

Fr = ms - bl;

figure
scatter(Dir,Fr,25,'b','filled')
xticks([0 45 90 135 180 225 270 315])
xticklabels({'0','45','90','135','180','225','270','315'})
xlim([-10 325])
% ylim([min(Fr)-1 max(Fr)+1])
ylim([min(Fr)-1 7])
title(type)
xlabel('Preferred Direction')
ylabel({'\deltaFr = Fr_{MS} - Fr_{NMS}'})

