%%
clear all; clc
cd E:\MT_MST\Microstim\PSTH\zscore
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

% MST.ri_MS_T = x(3).riMS;
% MST.ri_NMS_T = x(3).riNMS;
% MST.ri_MS_S = x(4).riMS;
% MST.ri_NMS_S = x(4).riNMS;

% pNMS = [0.5 theta(i,2) theta(i,6)  min(yNMS(:,i))];
% x(z).pHatMS(i,:) = lsqcurvefit(@vonMises,pMS,dirs,yMS(:,i)',[0 0 0 0],[2*pi 2*pi theta(i,5)+min(yMS(:,i)) min(yMS(:,i))]);
%% Tuning curves: MT & MST : Translation & Spiral : No-MicroStim & MicroStim
load(['E:\MT_MST\Microstim\PSTH\' 'ctl.mat'])
names = {'ytu310a','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a'};
haha = [5;8;nan;nan;1;4;5];
info = struct;
for i = 1:size(names,2) 
    info(i).MT = find(ctl_MT(:,2)==i);
    info(i).MST = find(ctl_MST(:,2)==i);
    info(i).MTcell = load(['E:\MT_MST\Microstim\Cell_Lib\Candidate Cells\','CLib_',names{1,i}(1:end-1),'.mat']);
    info(i).MSTcell = load(['E:\MT_MST\Microstim\Cell_Lib\MS Cells\','CLib_',names{1,i}(1:end-1),'.mat']);
    count(i,1) = size(info(i).MT,1);
    count(i,2) = size(info(i).MST,1);
    count(i,3) = count(i,1) + count(i,2);
end
for i = 1:size(names,2)
    figure
    p = ceil(sqrt(count(i,3)));
    ct = 0;
    for j = 1:count(i,2)
        ct = ct + 1;
        subplot(p,p,ct)
        if j == haha(i) && ~isnan(haha(i))
            plot(rad2deg(dirs),MST.ri_NMS_T(:,info(i).MST(j)),'Color','b','LineWidth',3,'LineStyle','--'); hold on
            plot(rad2deg(dirs),MST.ri_MS_T(:,info(i).MST(j)),'Color','b','LineWidth',3)
            plot(rad2deg(dirs),MST.ri_NMS_S(:,info(i).MST(j)),'Color','r','LineWidth',3,'LineStyle','--'); hold on
            plot(rad2deg(dirs),MST.ri_MS_S(:,info(i).MST(j)),'Color','r','LineWidth',3)
            xlim([0 315])
        else
            plot(rad2deg(dirs),MST.ri_NMS_T(:,info(i).MST(j)),'Color','b','LineWidth',1,'LineStyle','--'); hold on
            plot(rad2deg(dirs),MST.ri_MS_T(:,info(i).MST(j)),'Color','b','LineWidth',1)
            plot(rad2deg(dirs),MST.ri_NMS_S(:,info(i).MST(j)),'Color','r','LineWidth',1,'LineStyle','--'); hold on
            plot(rad2deg(dirs),MST.ri_MS_S(:,info(i).MST(j)),'Color','r','LineWidth',1)
            xlim([0 315])
        end
        title(sprintf('MST-%d',info(i).MSTcell.CLib(j)))
    end
%     legend({'Trans--> NMStim','Trans--> MStim','Spiral--> NMStim','Spiral--> MStim'})
    for z = 1:count(i,1)
        ct = ct + 1;
        subplot(p,p,ct)
        plot(rad2deg(dirs),MT.ri_NMS_T(:,info(i).MT(z)),'Color','m','LineWidth',1,'LineStyle','--'); hold on
        plot(rad2deg(dirs),MT.ri_MS_T(:,info(i).MT(z)),'Color','m','LineWidth',1)
        plot(rad2deg(dirs),MT.ri_NMS_S(:,info(i).MT(z)),'Color','k','LineWidth',1,'LineStyle','--'); hold on
        plot(rad2deg(dirs),MT.ri_MS_S(:,info(i).MT(z)),'Color','k','LineWidth',1)
        xlim([0 315])
        title(sprintf('MT-%d',info(i).MTcell.CLib(z)))
    end
%     legend({'Trans--> NMStim','Trans--> MStim','Spiral--> NMStim','Spiral--> MStim'})
    xlabel 'Direction'; ylabel 'Fr (spk/sec)';
    suptitle(sprintf('Session %d',i))
end
%% Tuning curves: MT & MST : Translation & Spiral : No-MicroStim & MicroStim
load(['E:\MT_MST\Microstim\PSTH\' 'ctl.mat'])
MT = load(['E:\MT_MST\Microstim\PSTH\' 'von_MT.mat']);
MST = load(['E:\MT_MST\Microstim\PSTH\' 'von_MST.mat']);
names = {'ytu310a','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a'};
haha = [5;8;nan;nan;1;4;5];
dirs = linspace(-pi,pi,8+1);
info = struct;
for i = 1:size(names,2) 
    info(i).MT = find(ctl_MT(:,2)==i);
    info(i).MST = find(ctl_MST(:,2)==i);
    info(i).MTcell = load(['E:\MT_MST\Microstim\Cell_Lib\Candidate Cells\','CLib_',names{1,i}(1:end-1),'.mat']);
    info(i).MSTcell = load(['E:\MT_MST\Microstim\Cell_Lib\MS Cells\','CLib_',names{1,i}(1:end-1),'.mat']);
    count(i,1) = size(info(i).MT,1);
    count(i,2) = size(info(i).MST,1);
    count(i,3) = count(i,1) + count(i,2);
end
for i = 1:size(names,2)
    figure
    p = ceil(sqrt(count(i,3)));
    ct = 0;
    for j = 1:count(i,2)
        ct = ct + 1;
        subplot(p,p,ct)
        if j == haha(i) && ~isnan(haha(i))
            plot(rad2deg(dirs),MST.ri_NMS_T(:,info(i).MST(j)),'Color','b','LineWidth',3,'LineStyle','--'); hold on
            plot(rad2deg(dirs),MST.ri_MS_T(:,info(i).MST(j)),'Color','b','LineWidth',3)
            plot(rad2deg(dirs),MST.ri_NMS_S(:,info(i).MST(j)),'Color','r','LineWidth',3,'LineStyle','--'); hold on
            plot(rad2deg(dirs),MST.ri_MS_S(:,info(i).MST(j)),'Color','r','LineWidth',3)
            xlim([-185 185])
        else
            plot(rad2deg(dirs),MST.ri_NMS_T(:,info(i).MST(j)),'Color','b','LineWidth',1,'LineStyle','--'); hold on
            plot(rad2deg(dirs),MST.ri_MS_T(:,info(i).MST(j)),'Color','b','LineWidth',1)
            plot(rad2deg(dirs),MST.ri_NMS_S(:,info(i).MST(j)),'Color','r','LineWidth',1,'LineStyle','--'); hold on
            plot(rad2deg(dirs),MST.ri_MS_S(:,info(i).MST(j)),'Color','r','LineWidth',1)
            xlim([-185 185])
        end
        title(sprintf('MST-%d',info(i).MSTcell.CLib(j)))
    end
    legend({'Trans--> NMStim','Trans--> MStim','Spiral--> NMStim','Spiral--> MStim'})
    for z = 1:count(i,1)
        ct = ct + 1;
        subplot(p,p,ct)
        plot(rad2deg(dirs),MT.ri_NMS_T(:,info(i).MT(z)),'Color','m','LineWidth',1,'LineStyle','--'); hold on
        plot(rad2deg(dirs),MT.ri_MS_T(:,info(i).MT(z)),'Color','m','LineWidth',1)
        plot(rad2deg(dirs),MT.ri_NMS_S(:,info(i).MT(z)),'Color','k','LineWidth',1,'LineStyle','--'); hold on
        plot(rad2deg(dirs),MT.ri_MS_S(:,info(i).MT(z)),'Color','k','LineWidth',1)
        xlim([-185 185])
        title(sprintf('MT-%d',info(i).MTcell.CLib(z)))
    end
    legend({'Trans--> NMStim','Trans--> MStim','Spiral--> NMStim','Spiral--> MStim'})
    xlabel 'Direction'; ylabel 'Fr (spk/sec)';
    suptitle(sprintf('Session %d',i))
end