clear all; clc
cd E:\MT_MST\Microstim\PSTH\MS\
type = 'trans';
option = 'MS';
load(sprintf('S%s%s.mat',type,option))
S = info;
clear info
load(sprintf('N%s%s.mat',type,option))
N = info;
names = {'ytu310a','ytu323a','ytu333a'};
% S(2:3) = []; N(2:3) = [];
%% Vector average: Preferred Direction
spd = horzcat(S.df); % Micro Stim
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
save(['E:\MT_MST\Microstim\PSTH\MS\',sprintf('theta_%s_%s.mat',type(1),option)],'theta')
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
        save(sprintf('trMST_%s_%s.mat',type(1),option),'trMST','tr')        
end
%% PD difference between MST(NoMS) & MT(NoMS & MS): Part II
clear all; clc
load(['E:\MT_MST\Microstim\PSTH\MS\' 'trMST_t_MS.mat']); trMSTt = trMST; clear trMST
load(['E:\MT_MST\Microstim\PSTH\MS\' 'trMST_s_MS.mat']); trMSTs = trMST; clear trMST

load(['E:\MT_MST\Microstim\PSTH\MS\' 'theta_t_.mat']); trMTt = theta; clear theta
load(['E:\MT_MST\Microstim\PSTH\MS\' 'theta_s_.mat']); trMTs = theta; clear theta

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
save('E:\MT_MST\Microstim\PSTH\MS\pdDiff.mat','pdDiff')
save('E:\MT_MST\Microstim\PSTH\MS\pdD.mat','pdD')
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
cd E:\MT_MST\Microstim\PSTH\MS\
option = 'MS';
load(sprintf('theta_t_%s.mat',option)); theta_T = theta; clear theta
load(sprintf('theta_s_%s.mat',option)); theta_S = theta; clear theta
figure
subplot(2,2,1)
polarhistogram(theta_T(:,3),20,'FaceColor','b','EdgeColor','g','FaceAlpha',0.7,'EdgeAlpha',.5)
title 'P^{Trans}_{Dir}(S-N)';
subplot(2,2,2)
polarhistogram(theta_S(:,3),20,'FaceColor','k','EdgeColor','r','FaceAlpha',0.7,'EdgeAlpha',.5)
title 'P^{Spiral}_{Dir}(S-N)';
subplot(2,2,3)
scatter(theta_T(:,5),theta_T(:,6),'MarkerEdgeColor','g','MarkerFaceColor','b'); refline(1,0)
xlim([0 max(theta_T(:,5))+10])
ylim([0 max(theta_T(:,6))+10])
title 'R^{Trans}_{peak}'
xlabel 'MicroStim(Hz)' 
ylabel 'No MicroStim(Hz)'
subplot(2,2,4)
scatter(theta_S(:,5),theta_S(:,6),'MarkerEdgeColor','r','MarkerFaceColor','k'); refline(1,0)
xlim([0 max(theta_S(:,5))+10])
ylim([0 max(theta_S(:,6))+10])
title 'R^{Spiral}_{peak}'
xlabel 'MicroStim(Hz)' 
ylabel 'No MicroStim(Hz)'
%% Normalized Response as a Function of Position - based on max firing rate - RF Center of Mass - RF Shift
% Part I: Max firing rate btw 8 dir @ each 9 pos
% Part II: Center of Mass
% Part III: RF Center of Mass
% Part IV: RF Shift

% Max firing rate between 8 directions @ each 9 position
posInf = struct;
for i = 1:size(S,2)
    for j = 1:size(S(i).xr,2)
        posInf(i).fmaxposS(j,:) = max(S(i).xr(j).firing); % Max firing rate between 8 directions @ each 9 position
        posInf(i).fmaxposN(j,:) = max(N(i).xr(j).firing);
    end
end
fmaxS = vertcat(posInf.fmaxposS);
fmaxN = vertcat(posInf.fmaxposN);
fmaxS(badCell,:) = [];
fmaxN(badCell,:) = [];
frmax(:,1) = max(fmaxS');
frmax(:,2) = max(fmaxN');
%%%% Center of Mass %%%%
sepration = 17;
diam = 20;
x = [-17;0;17;-17;0;17;-17;0;17]; 
y = [17;17;17;0;0;0;-17;-17;-17];
for i = 1:size(fmaxS,1)
    [Vx(i,1) Vy(i,1)] = cmass(fmaxS(i,:)',x,y);
    [Vx(i,2) Vy(i,2)] = cmass(fmaxN(i,:)',x,y);
end
diffV(:,1) = Vx(:,1) - Vx(:,2); % x-axis
diffV(:,2) = Vy(:,1) - Vy(:,2); % y-axis
%%%% RF center of mass %%%%
figure
scatter(Vx(:,1),Vy(:,1),'FaceColor','b'); hold on
scatter(Vx(:,2),Vy(:,2),'FaceColor','r') 
line([0 0],[min(y) max(y)],'LineWidth',0.5,'Color','k')
line([min(x) max(x)],[0 0],'LineWidth',0.5,'Color','k')
xlim([min(x) max(x)])
ylim([min(y) max(y)])
title(sprintf('RF^{%s}_{cm}',type))
xlabel 'Position^o'
ylabel 'Position^o'
legend({'MicroStim','No-MicroStim'},'location','NorthWest')
%%%% RF shift %%%%
figure
for i = 1:size(diffV,1) 
    quiver(Vx(i,2),Vy(i,2),diffV(i,1),diffV(i,2),'LineWidth',1,'Color','k','MaxHeadSize',.6); hold on 
end
line([0 0],[min(y) max(y)],'LineWidth',0.1,'Color','g')
line([min(x) max(x)],[0 0],'LineWidth',0.1,'Color','g')
xlim([min(x) max(x)])
ylim([min(y) max(y)])
title(sprintf('RF^{%s}_{cm}',type))
xlabel 'Position^o'
ylabel 'Position^o'
legend({'Shift'},'location','NorthWest')
%%%% Net arrow calculation + Plot %%%%
figure
for i = 1:size(diffV,1) 
    quiver(0,0,diffV(i,1),diffV(i,2),'LineWidth',1,'Color','k','MaxHeadSize',.6); hold on 
end
xx = sum(diffV(:,1)); % x-axis
yy = sum(diffV(:,2)); % y-axis
quiver(0,0,xx,yy,'LineWidth',1.5,'Color','b','MaxHeadSize',.6); hold on
line([0 0],[min(y) max(y)],'LineWidth',0.1,'Color','g')
line([min(x) max(x)],[0 0],'LineWidth',0.1,'Color','g')
% xlim([-abs(xx) abs(xx)])
% ylim([-abs(yy) abs(yy)])
xlim([-5 5])
ylim([-5 5])
% ylim([-abs(yy) abs(yy)])
title(sprintf('RF^{%s}_{cm}',type))
xlabel 'Position^o'
ylabel 'Position^o'
return
%% Tuning curves: MT & MST : Translation & Spiral : No-MicroStim & MicroStim
clear all; clc
load('E:\MT_MST\Microstim\PSTH\MS\ctl_MT.mat')
load('E:\MT_MST\Microstim\PSTH\MS\MT.mat')
load('E:\MT_MST\Microstim\PSTH\MS\pdDiff.mat')
load('E:\MT_MST\Microstim\PSTH\MS\pdD.mat')
load('E:\MT_MST\Microstim\PSTH\MS\trMST_t_MS.mat','tr'); trt = tr; clear tr
load('E:\MT_MST\Microstim\PSTH\MS\trMST_s_MS.mat','tr'); trs = tr; clear tr

names = {'ytu310a','ytu323a','ytu333a'};
dirs = deg2rad(0:45:315);
info = struct;
for i = 1:size(names,2) 
    info(i).MT = find(ctl_MT(:,2)==i);
    info(i).PDdiffT = pdDiff(info(i).MT,2); % MST(NoMS) - MT(NoMS) for Translation
    info(i).PDdiffS = pdDiff(info(i).MT,4); % MST(NoMS) - MT(NoMS) for Spiral
    info(i).PDdT = pdD(info(i).MT,1); % MT(NoMS) - MT(MS) for Translation
    info(i).PDdS = pdD(info(i).MT,2); % MT(NoMS) - MT(MS) for Spiral
    info(i).MTcell = load(['E:\MT_MST\Microstim\Cell_Lib\Candidate Cells\','CLib_',names{1,i}(1:end-1),'.mat']);
    info(i).MSTcell = load(['E:\MT_MST\Microstim\Cell_Lib\MS Cells\MUA-cleaned\','CLib_',names{1,i}(1:end-1),'.mat']);
    count(i,1) = size(info(i).MT,1);
    count(i,2) = size(trt(i).firingS,2);
    count(i,3) = count(i,1) + count(i,2);
end
%%%% Plot section %%%%
for i = 1:size(names,2)
    figure
    p = ceil(sqrt(count(i,3)));
    ct = 0;
    for j = 1:count(i,2)
        ct = ct + 1;
        subplot(p,p,ct)
%         plot(rad2deg(dirs),trt(i).firingN(:,j),'Color','m','LineWidth',1,'LineStyle','--'); hold on
        plot(rad2deg(dirs),trs(i).firingN(:,j),'Color','m','LineWidth',1,'LineStyle','--'); hold on
        xlim([0 315])
        title(sprintf('MST-%d',info(i).MSTcell.CLib(j)))
    end
    for z = 1:count(i,1)
        ct = ct + 1;
        subplot(p,p,ct)
%         plot(rad2deg(dirs),MT.ri_NMS_T(:,info(i).MT(z)),'Color','b','LineWidth',1,'LineStyle','--'); hold on
%         plot(rad2deg(dirs),MT.ri_MS_T(:,info(i).MT(z)),'Color','b','LineWidth',2)
        plot(rad2deg(dirs),MT.ri_NMS_S(:,info(i).MT(z)),'Color','k','LineWidth',1,'LineStyle','--'); hold on
        plot(rad2deg(dirs),MT.ri_MS_S(:,info(i).MT(z)),'Color','k','LineWidth',2)
        xlim([0 315])
%         title(sprintf('{MT_{NoMS}- MT_{MS}} = %.f^o / {MST_{NoMS}- MT_{NoMS}} = %.f^o',info(i).PDdT(z,1),info(i).PDdiffT(z,1)))
        title(sprintf('{MT_{NoMS}- MT_{MS}} = %.f^o / {MST_{NoMS}- MT_{NoMS}} = %.f^o',info(i).PDdS(z,1),info(i).PDdiffS(z,1)))
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