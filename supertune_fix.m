clear all
close all
clc
name = 'ytu337c';
load(['E:\MT_MST\Microstim\Cell_Lib\MST_MU_MT_SU\MT\','CLib_',name(1:end-1),'.mat'])
chunum = CLib;
load(['E:\MT_MST\Plexon\RFiles\',name,'_TrialStructure.mat'])
Fs = 10000;
sv = 1;
nev = 0;
flag = 1; %1 or 3
% MS timing properties
FsMS = 200;
nPulse = 5;
PulseDur = floor((1/FsMS) * nPulse * 1000);
OnsetMS = 0.150; 
tstMS = OnsetMS * 1000;
tendMS = tstMS + PulseDur;
%
for ci = 1:length(chunum)
    counting = 1;
    ch = chunum(ci);
    unit = 1;
    if nev == 1
        [spikeMat,numspikes,stimLength,baseLineLength] = nev2NFile_superTune(name,ch,unit);
    elseif nev == 0
        [spikeMat,numspikes,stimLength,baseLineLength] = plx2NFile_superTune(name,ch,unit,flag);
    end
    angles = sort(unique(spikeMat(:,3)));
    S = angles(logical([1 0 1 0 1 0 1 0]));
    N = angles(~logical([1 0 1 0 1 0 1 0]));
    for i = 1:length(spikeMat)
        if any(spikeMat(i,3) == S)
            spikeMat(i,3) = 1;
        elseif any(spikeMat(i,3) == N)
            spikeMat(i,3) = 0;
        end
    end
    spikeMat(:,4) = 1; % Motion type
    spikeMat(:,7) = 1; % Position
    t1 = 0;
    t2 = mean(stimLength);
    time = floor((t2-t1)*Fs);
    time_bl = floor(mean(baseLineLength)*Fs);
    featureIdx = 3; %the 8 directions for all motion types
    xaxis = unique(spikeMat(:,featureIdx));
    numMotion = max(unique(spikeMat(:,4)));
    
    vec = zeros(length(unique(spikeMat(:,featureIdx))),1);
    startstim = 0;
    endstim = mean(stimLength) * 1000;
    firing = [];
    firing_bl = [];
    motionType = numMotion;
    spikes = spikeMat(spikeMat(:,1)>startstim & spikeMat(:,4)==motionType,:);
    spikes(0<spikes(:,1) & spikes(:,1)<tendMS,:) = [];
    trialsPerFeature = floor(length(unique(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==motionType,2))) / length(unique(spikeMat(:,3))));
    if counting == 1
        spktrain = zeros(time,trialsPerFeature);
        spktrain_bl = zeros(time_bl,trialsPerFeature);
    end
    for j = 1:length(xaxis)
        vec = length(spikes(spikes(:,featureIdx)==xaxis(j),1));
        vec_rep_sep = spikes(spikes(:,featureIdx)==xaxis(j),:);
        spikes_bl = spikeMat(spikeMat(:,1)<=0 & spikeMat(:,1)>=(0-time_bl)*1000/Fs & spikeMat(:,4)==motionType & spikeMat(:,featureIdx)==xaxis(j),:);
        trialsPerThisFeature = floor(length(unique(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==motionType & spikeMat(:,featureIdx)==xaxis(j),2))));
        firing = vec / (trialsPerThisFeature * ((endstim - tendMS) - startstim) / 1000);
        firing_bl = length(spikes_bl(:,1))/(trialsPerThisFeature * (mean(baseLineLength)));
        [RFcenter,I] = max(firing(:));
        [i1,i2,RFcenterIdx,i3,i4] = ind2sub(size(firing),I);
        fr(ci,1,j) = firing;
        fr(ci,2,j) = firing_bl;
        switch xaxis(j)
            case 0
                type = 'N';
            case 1
                type = 'S';
        end
        if sv == 1
            save(['E:\MT_MST\SuperTuneFiringMatrix\fixation\',name,num2str(ch),num2str(unit),type,'firingMat.mat'],'firing')
            save(['E:\MT_MST\SuperTuneFiringMatrix\fixation\',name,num2str(ch),num2str(unit),type,'firingMat_bl.mat'],'firing_bl')
        end
    end
end
save(['E:\MT_MST\SuperTuneFiringMatrix\fixation\',name,'fr.mat'],'fr')
return
%%
clear all; clc
names = {'ytu310a','ytu316a','ytu321a','ytu329a','ytu333a','ytu335a','ytu337a'};
frn = []; frs = [];
for i = 1:length(names)
    fn = names{i};
    a = struct2array(load(['E:\MT_MST\SuperTuneFiringMatrix\fixation\',fn(1:end-1),'cfr.mat']));
    frn = [frn;a(:,:,1)]; %N
    frs = [frs;a(:,:,2)]; %S
end
%%
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
%%
% for o = 1:trialsPerThisFeature
%     trialIdx = unique(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==motionType & spikeMat(:,featureIdx)==xaxis(j),2));
%     spktrain1 = zeros(time,1);
%     indices = round(vec_rep_sep(vec_rep_sep(:,2)==trialIdx(o),1)*(Fs/1000));
%     indices = indices(indices<=time);
%     spktrain1(indices) = 1;
%     spktrain(:,o) = spktrain1;
%     spktrain2 = zeros(time_bl,1);
%     if trialIdx(o)> length(baseLineLength)
%         indices_bl = round((spikes_bl(spikes_bl(:,2)==trialIdx(o),1)+(mean(baseLineLength)*1000))*(Fs/1000));
%     else
%         indices_bl = round((spikes_bl(spikes_bl(:,2)==trialIdx(o),1)+(baseLineLength(trialIdx(o))*1000))*(Fs/1000));
%     end
%     indices_bl = indices_bl((indices_bl>0)&(indices_bl<=time_bl));
%     spktrain2(indices_bl) = 1;
%     spktrain_bl(:,o) = spktrain2;
% end
% save(['E:\MT_MST\SuperTuneSpkTrains\fixation\',name,num2str(ch),num2str(unit),'spktrain_bl.mat'],'spktrain_bl','Fs')
% save(['E:\MT_MST\SuperTuneSpkTrains\fixation\',name,num2str(ch),num2str(unit),'spktrain.mat'],'spktrain','Fs','RFcenterIdx')