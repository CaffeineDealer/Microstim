clear all
close all
clc
Fs = 10000;
sv = 1;
flag = 1; %1 or 3
name = 'ytu337b';
FsMS = 200;
nPulse = 5;
PulseDur = floor((1/FsMS) * nPulse * 1000);
OnsetMS = 0.150; 
tstMS = OnsetMS * 1000;
tendMS = tstMS + PulseDur;
% load(['E:\MT_MST\Microstim\Cell_Lib\Candidate Cells\','CLib_',name(1:end-1),'.mat'])
CLib = 1:64';
chunum = CLib;
load(['E:\MT_MST\Plexon\RFiles\',name,'_TrialStructure.mat'])
rows = file.taskDialogValues.superTuneRows;
columns = file.taskDialogValues.superTuneColumns;
ci = 1; counting = 1; ch = 1; unit = 1;
% load('E:\MT_MST\Plexon\RFiles\nevlib.mat')
%%
for ci = 1:length(chunum)
    counting = 1;
    ch = chunum(ci);
    unit = 1;
    [spikeMat,numspikes,stimLength,baseLineLength] = plx2NFile_superTune(name,ch,unit,flag);
    t1 = 0;
    t2 = mean(stimLength);
    time = floor((t2-t1)*Fs);
    time_bl = floor(mean(baseLineLength)*Fs);
    featureIdx = 3; %the 8 directions for all motion types
    xaxis = unique(spikeMat(:,featureIdx));
    numMotion = max(unique(spikeMat(:,4)));
    sizes = unique(spikeMat(:,8));
    count = 1:length(xaxis);
    vec = zeros(length(unique(spikeMat(:,featureIdx))),1);
    gridIndeces = unique(spikeMat(:,7));
    coherences = unique(spikeMat(:,9));
    numTrials = length(spikeMat(spikeMat(:,1)==-1000,1));
    startstim = 0;
    endstim = mean(stimLength) * 1000;
    firing = zeros(length(xaxis),numMotion,rows*columns,length(sizes),length(coherences));
    firing_bl = firing;
    responses1 = zeros(length(xaxis),numMotion,rows,columns);
    for m = 1:length(coherences)
        coh = coherences(m);
        for j = 1:length(sizes)
            siz = sizes(j);
            for l = 1:numMotion
                motionType = l;
                for k = 1:max(gridIndeces)
                    gridIndex = k;
                    spikes = spikeMat(spikeMat(:,1)>startstim & spikeMat(:,4)==motionType & spikeMat(:,7)==gridIndex & spikeMat(:,8)==siz & spikeMat(:,9)==coh,:);
%                     spikes(tstMS<spikes(:,1) & spikes(:,1)<tendMS,:) = [];
                    spikes(0<spikes(:,1) & spikes(:,1)<tendMS,:) = [];
                    spikes(spikes(:,1)>300,:) = [];
%                     spikes(tstMS<spikes(:,1),:) = [];
                    trialsPerFeature = floor(length(unique(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==motionType & spikeMat(:,7)==gridIndex & spikeMat(:,8)==siz & spikeMat(:,9)==coh,2)))...
                        /length(unique(spikeMat(:,3))));
                    if counting == 1
                        spktrain = zeros(time,length(xaxis),numMotion,rows*columns,trialsPerFeature,length(sizes),length(coherences));
                        spktrain_bl = zeros(time_bl,length(xaxis),numMotion,rows*columns,trialsPerFeature,length(sizes),length(coherences));
                        counting = 2;
                    end
                    for i = count
                        vec(i) = length(spikes(spikes(:,featureIdx)==xaxis(i),1));
                        vec_rep_sep = spikes(spikes(:,featureIdx)==xaxis(i),:);
                        spikes_bl = spikeMat(spikeMat(:,1)<=0 & spikeMat(:,1)>=(0-time_bl)*1000/Fs & spikeMat(:,4)==motionType &...
                            spikeMat(:,7)==gridIndex & spikeMat(:,featureIdx)==xaxis(i) & spikeMat(:,8)==siz & spikeMat(:,9)==coh,:);
                        trialsPerThisFeature = floor(length(unique(spikeMat(spikeMat(:,1)==-1000 & ...
                            spikeMat(:,4)==motionType & spikeMat(:,7)==gridIndex  & spikeMat(:,featureIdx)==xaxis(i)&...
                            spikeMat(:,8)==siz & spikeMat(:,9)==coh,2))));
                        firing_bl(i,l,k,j,m) = length(spikes_bl(:,1))/(trialsPerThisFeature*(mean(baseLineLength)));
                        if trialsPerThisFeature > trialsPerFeature
                            trialsPerThisFeature = trialsPerFeature;
                        end
                        for o = 1:trialsPerThisFeature
                            trialIdx = unique(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==motionType & spikeMat(:,7)==gridIndex & spikeMat(:,featureIdx)==xaxis(i)& spikeMat(:,8)==siz & spikeMat(:,9)==coh,2));
                            spktrain1 = zeros(time,1);
                            indices = round(vec_rep_sep(vec_rep_sep(:,2)==trialIdx(o),1)*(Fs/1000));
                            indices = indices(indices<=time);
                            spktrain1(indices) = 1;
                            spktrain(:,i,l,k,o,j,m) = spktrain1;
                            spktrain2 = zeros(time_bl,1);
                            if trialIdx(o)> length(baseLineLength)
                                indices_bl = round((spikes_bl(spikes_bl(:,2)==trialIdx(o),1)+(mean(baseLineLength)*1000))*(Fs/1000));
                            else
                                indices_bl = round((spikes_bl(spikes_bl(:,2)==trialIdx(o),1)+(baseLineLength(trialIdx(o))*1000))*(Fs/1000));
                            end
                            indices_bl = indices_bl((indices_bl>0)&(indices_bl<=time_bl));
                            spktrain2(indices_bl) = 1;
                            spktrain_bl(:,i,l,k,o,j,m) = spktrain2;
                        end
%                         firing(i,l,k,j,m) = vec(i) / (trialsPerThisFeature * ((endstim - tendMS) - startstim) / 1000);
                        firing(i,l,k,j,m) = vec(i) / (trialsPerThisFeature * ((300 - tendMS) - startstim) / 1000);
                    end
                end
                [RFcenter,I] = max(firing(:));
                [i1,i2,RFcenterIdx,i3,i4] = ind2sub(size(firing),I);
            end
        end
    end
    if sv == 1
        save(['E:\MT_MST\SuperTuneSpkTrains\ms_rmvd_150300\',name,num2str(ch),num2str(unit),'spktrain_bl.mat'],'spktrain_bl','Fs')
        save(['E:\MT_MST\SuperTuneSpkTrains\ms_rmvd_150300\',name,num2str(ch),num2str(unit),'spktrain.mat'],'spktrain','Fs','RFcenterIdx')
        save(['E:\MT_MST\SuperTuneFiringMatrix\ms_rmvd_150300\',name,num2str(ch),num2str(unit),'firingMat.mat'],'firing')
        save(['E:\MT_MST\SuperTuneFiringMatrix\ms_rmvd_150300\',name,num2str(ch),num2str(unit),'firingMat_bl.mat'],'firing_bl')
    end
    clearvars -except plt sv v3categ flag nev DprimeV3 DprimeMT RFmapping sep sizeTuning firingAve channel firingGp chunum name Fs microstim condition rows columns microstimParadigm tstMS tendMS
end