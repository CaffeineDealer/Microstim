clear all; clc
cd E:\MT_MST\Microstim\PSTH\msr\
type = 'spiral';
option = '';
load(sprintf('S%s%s.mat',type,option))
S = info;
clear info
load(sprintf('N%s%s.mat',type,option))
N = info;
switch option
    case ''
        names = {'ytu310a','ytu312b','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a'};
%         S(4) = []; N(4) = [];
%         S(6) = []; N(6) = [];
    case 'MS'
        names = {'ytu310a','ytu312b','ytu316a','ytu323a','ytu333a','ytu335a','ytu337a'};
%         S(7) = []; N(7) = [];
end
%% Vector average: Preferred Direction
% spd = horzcat(S.df); % Micro Stim
n = 1;
for i = 1:size(N,2)
    for j = 1:size(S(i).xr,2)
        spd(:,n) = S(i).xr(j).firing(:,N(i).pt(j));
        n = n + 1;
    end
end
npd = horzcat(N.df); % No Micro Stim
%%%% Form a library for cell ID and File ID + change/no change in PD %%%%
ctl = vertcat(S.ch); % Ch ID = ctl(:,1)
fid = vertcat(S.fid);
ct = 0;
for i = 1:size(fid,1)
    if i == 1 
        ctl(1:fid(i,1),2) = i; % File ID = ctl(:,2)
    else
        ctl(ct+1:ct+fid(i,1),2) = i;
    end
     ct = ct + fid(i,1);
end
%%%% Vector Average %%%%
theta = nan(size(spd,2),6);
r = nan(size(spd,2),1);
x = nan(8,size(spd,2));
for i = 1:size(spd,2)
    if ~isnan(spd(1,i)) & ~isnan(npd(1,i))
        [theta(i,1),~,r(i,1)] = polarplot(spd(:,i)); % Micro Stim
        [theta(i,2),~,r(i,1)] = polarplot(npd(:,i)); % No Micro Stim
        theta(i,3) = abs(angdiff(theta(i,1),theta(i,2))); % Degree difference
        theta(i,4) = i; % cell position
        theta(i,5) = max(spd(:,i)); % Max Firing Rate for Micro Stim
        theta(i,6) = max(npd(:,i)); % Max Firing Rate for No Micro Stim
        theta(i,7) = ctl(i,2); % cell session
        x(:,i) = spd(:,i);
        if theta(i,3) > deg2rad(20) % Define a change/no change criteria
            ctl(i,3) = 1; % Cells with change in PD
            ctl(i,6) = rad2deg(theta(i,3));
        else
            ctl(i,3) = 0; % Cells with no change in PD
            ctl(i,6) = rad2deg(theta(i,3));
        end
    else
        ctl(i,3) = nan; % Bad Cells
    end
end
badCell = find(isnan(theta(:,1)));
ctl(:,4) = vertcat(S.pt); % Micro stim spatial position
ctl(:,5) = vertcat(N.pt); % No Micro stim spatial position
theta(isnan(theta(:,1)),:) = [];
r(isnan(r(:,1))) = [];
[rv, pv] = corrcoef(theta(:,1),theta(:,2));
save(['E:\MT_MST\Microstim\PSTH\msr\',sprintf('theta_%s_%s.mat',type(1),option)],'theta')
save(['E:\MT_MST\Microstim\PSTH\msr\',sprintf('mt_%s%s.mat',type,option)],'theta','spd','npd')
RF(S,N,badCell,type,theta)
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
a = zeros(77,14);
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
% pos = vertcat(N.pt);
% dir = vertcat(N.Dir); dir(:,2) = [];
n = 1;
for i = 1:size(N,2)
    for j = 1:size(S(i).xr,2)
        spd(:,n) = S(i).xr(j).firing(:,N(i).pt(j));
        n = n + 1;
    end
end
[frm frpos] = max(npd);
for i = 1:size(frpos,2)
    frm(2,i) = spd(frpos(i),i);
end
frm(3,:) = max(spd);
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
        scatter(frm(1,:),frm(2,:),'MarkerEdgeColor','g','MarkerFaceColor','b'); refline(1,0)
    case 'spiral'
        cmap(frm(4,:)==1,1) = 1;
        scatter(frm(1,:),frm(2,:),'MarkerEdgeColor','r','MarkerFaceColor','k'); refline(1,0)
end
% scatter(frm(1,:),frm(2,:),25,cmap,'filled'); refline(1,0)
title(sprintf('Resp^{%s}_{MT}',type))
xlim([0 100]); xlabel 'Fr^{NMS} (spk/sec)' % max(frm(1,:))
ylim([0 90]); ylabel 'Fr^{MS} (spk/sec)' % max(frm(2,:))
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
%% PD difference between MST(NoMS) & MT(NoMS & MS): Part II
clear all; clc
load(['E:\MT_MST\Microstim\PSTH\msr\mt-mst\' 'trMST_t_MS.mat']); trMSTt = trMST; clear trMST
load(['E:\MT_MST\Microstim\PSTH\msr\mt-mst\' 'trMST_s_MS.mat']); trMSTs = trMST; clear trMST

load(['E:\MT_MST\Microstim\PSTH\msr\mt-mst\' 'theta_t_.mat']); trMTt = theta; clear theta
load(['E:\MT_MST\Microstim\PSTH\msr\mt-mst\' 'theta_s_.mat']); trMTs = theta; clear theta

pdDiff = zeros(size(trMTt,1),4);
for i = 1:max(trMTt(:,7))
    mtidx = find(trMTt(:,7)==i);
    mstT = ones(size(mtidx,1),1) * trMSTt(i,2);
    mstS = ones(size(mtidx,1),1) * trMSTs(i,2);
    pdDiff(mtidx,1) = abs(angdiff(mstT,trMTt(mtidx,1))); % MST Trans NoMS vs. MT Trans MS
    pdDiff(mtidx,2) = abs(angdiff(mstT,trMTt(mtidx,2))); % MST Trans NoMS vs. MT Trans NoMS
    pdDiff(mtidx,3) = abs(angdiff(mstS,trMTs(mtidx,1))); % MST Spiral NoMS vs. MT Spiral MS
    pdDiff(mtidx,4) = abs(angdiff(mstS,trMTs(mtidx,2))); % MST Spiral NoMS vs. MT Spiral NoMS
end
pdDiff = rad2deg(pdDiff);
pdD(:,1) = rad2deg(trMTt(:,3));
pdD(:,2) = rad2deg(trMTs(:,3));
save('E:\MT_MST\Microstim\PSTH\msr\mt-mst\pdDiff.mat','pdDiff')
save('E:\MT_MST\Microstim\PSTH\msr\mt-mst\pdD.mat','pdD')
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
cd E:\MT_MST\Microstim\PSTH\msr\
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
%% Normalized Response as a Function of Position - based on max firing rate - RF Center of Mass - RF Shift
% Part I: Max firing rate btw 8 dir @ each 9 pos
% Part II: Center of Mass
% Part III: RF Center of Mass
% Part IV: RF Shift

% Max firing rate between 8 directions @ each 9 position
RF(S,N,badCell,type,theta)
%% Tuning curves: MT & MST : Translation & Spiral : No-MicroStim & MicroStim
clear all; clc
load('E:\MT_MST\Microstim\PSTH\msr\mt-mst\ctl_MT.mat')
load('E:\MT_MST\Microstim\PSTH\msr\mt-mst\MT.mat')
load('E:\MT_MST\Microstim\PSTH\msr\mt-mst\pdDiff.mat')
load('E:\MT_MST\Microstim\PSTH\msr\mt-mst\pdD.mat')
load('E:\MT_MST\Microstim\PSTH\msr\mt-mst\trMST_t_MS.mat','tr'); trt = tr; clear tr
load('E:\MT_MST\Microstim\PSTH\msr\mt-mst\trMST_s_MS.mat','tr'); trs = tr; clear tr
type = 'trans';
names = {'ytu310a','ytu312b','ytu316a','ytu323a','ytu333a','ytu335a'};
dirs = deg2rad(0:45:315);
info = struct;
for i = 1:size(names,2) 
    info(i).MT = find(ctl_MT(:,2)==i);
    info(i).PDdiffT = pdDiff(info(i).MT,2); % MST(NoMS) - MT(NoMS) for Translation
    info(i).PDdiffS = pdDiff(info(i).MT,4); % MST(NoMS) - MT(NoMS) for Spiral
    info(i).PDdT = pdD(info(i).MT,1); % MT(NoMS) - MT(MS) for Translation
    info(i).PDdS = pdD(info(i).MT,2); % MT(NoMS) - MT(MS) for Spiral
    info(i).MTcell = load(['E:\MT_MST\Microstim\Cell_Lib\MS_removed\Candidate Cells\','CLib_',names{1,i}(1:end-1),'.mat']);
    info(i).MSTcell = load(['E:\MT_MST\Microstim\Cell_Lib\MS_removed\MS Cells\','CLib_',names{1,i}(1:end-1),'.mat']);
    count(i,1) = size(info(i).MT,1);
    count(i,2) = size(trt(i).firingS,2);
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
cd E:\MT_MST\Microstim\PSTH\msr\mt-mst
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
save(['E:\MT_MST\Microstim\PSTH\msr\mt-mst\MT.mat'],'MT')