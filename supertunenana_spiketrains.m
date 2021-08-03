clear all
close all
clc
Fs = 10000;
plt = 1; % Disable tuning plots
sv = 0;
microstim = 1;
nev = 1;
microstimParadigm = 0;
flag = 1; %1 or 3
RFmapping = 0;
condition = 'N';
name = 'slu066a';
FsMS = 200; 
nPulse = 5;
PulseDur = (1/FsMS) * nPulse * Fs; 
OnsetMS = 0.150;
tstMS = OnsetMS * Fs;
tendMS = tstMS + PulseDur;
unit = 1;
load(['E:\MT_MST\Synchrony things\chunum\',name(1:end-1),'Chunum.mat'])
% load(['E:\MT_MST\Microstim\Cell_Lib\MS Cells\new\','CLib_',name(1:end-1),'.mat'])
% load(['E:\MT_MST\Microstim\Cell_Lib\Candidate Cells\','CLib_',name(1:end-1),'.mat'])
% load(['E:\MT_MST\Microstim\Cell_Lib\MS Cells\','CLib_',name(1:end-1),'.mat'])
% chunum = CLib;
load(['E:\MT_MST\Plexon\RFiles\',name,'_TrialStructure.mat'])
rows = file.taskDialogValues.superTuneRows;
columns = file.taskDialogValues.superTuneColumns;
sep = file.taskDialogValues.superTuneSeperation;
ci = 1; counting = 1; ch = 1; unit = 1;
chunum = [27;28]
%%
for ci = 1:length(chunum)
    counting = 1;
    ch = chunum(ci);
    unit = 1;
    if microstim == 1
        load(['E:\MT_MST\Plexon\RFiles\',name(1:end-1),condition,num2str(ch),num2str(unit),'N.mat']);
%         load(['E:\MT_MST\Plexon\RFiles\sl1\',name(1:end-1),condition,num2str(ch),num2str(unit),'N.mat']);
        if strcmp(condition,'N')
            spikeMat = spikeMatnostim;
        elseif strcmp(condition,'S')
            spikeMat = spikeMatstim;
        elseif strcmp(condition,'Nf')
            spikeMat = spikeMatnostim;
        elseif strcmp(condition,'Sf')
            spikeMat = spikeMatstim;
        end
        stimInfo = load(['E:\MT_MST\MonkeyLab\HyperFlows\',name,'\',name,num2str(ch),num2str(unit),'N.mat']);
        stimLength = stimInfo.stimLength;
        baseLineLength = stimInfo.baseLineLength;
    elseif microstim == 0
        if nev == 1
            [spikeMat,numspikes,stimLength,baseLineLength] = nev2NFile_superTune(name,ch,unit);
        elseif nev == 0
            [spikeMat,numspikes,stimLength,baseLineLength] = plx2NFile_superTune(name,ch,unit,flag);
        end
    end
    t1 = 0;
    t2 = mean(stimLength);
    time = floor((t2-t1)*Fs);
    time_bl = floor(mean(baseLineLength)*Fs);
    %spikemat C2 trial#, C3 dir, C4 motion, C5 speed, spike times are aligned from stim on. 0 means stim on
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
                figure
                a = zeros(max(gridIndeces),4);
                for k = 1:max(gridIndeces)
                    gridIndex = k;
                    spikes = spikeMat(spikeMat(:,1)>startstim & spikeMat(:,4)==motionType & ...
                        spikeMat(:,7)==gridIndex & spikeMat(:,8)==siz & spikeMat(:,9)==coh,:);
                    trialsPerFeature = floor(length(unique(spikeMat(spikeMat(:,1)==-1000 & ...
                        spikeMat(:,4)==motionType & spikeMat(:,7)==gridIndex & spikeMat(:,8)==siz & spikeMat(:,9)==coh,2)))...
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
                        firing_bl(i,l,k,j,m) = length(spikes_bl(:,1))/...
                            (trialsPerThisFeature*(mean(baseLineLength)));
                        if trialsPerThisFeature > trialsPerFeature
                            trialsPerThisFeature = trialsPerFeature;
                        end
                        for o = 1:trialsPerThisFeature
                            trialIdx = unique(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==motionType &...
                                spikeMat(:,7)==gridIndex & spikeMat(:,featureIdx)==xaxis(i)& spikeMat(:,8)==siz & spikeMat(:,9)==coh,2));
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
                        firing(i,l,k,j,m) = vec(i)/(trialsPerThisFeature*(endstim-startstim)/1000);
                        spktrainAdj = spktrain;
                        spktrainAdj_bl = spktrain_bl;
                        spktrainAdj(tstMS:tendMS,:,:,:,:,:,:) = [];
                        spktrainAdj_bl(tstMS:tendMS,:,:,:,:,:,:) = [];
                        
                    end
                    nana(k) = subplot(rows,columns,k); %was 3 l k
                    plot(nana(k),xaxis*180/pi,squeeze(firing(:,l,k,j,m)));
                    hold on
%                     plot(xaxis*180/pi,firing_bl(:,l,k,j,m),'r')
                    f_bl_mean(1:size(xaxis,1)) = mean(firing_bl(:,l,k,j,m)); % Hacked by Yavar
                    plot(xaxis*180/pi,f_bl_mean,'r')
                    a(k,:) = axis();
                end
%                 f_bl_mean(1:size(xaxis,1)) = mean(firing_bl(:));
%                 for l = 1%:numMotion
%                     for k = 1:max(gridIndeces)
%                         nana(k) = subplot(rows,columns,k);
%                         plot(nana(k),xaxis*180/pi,squeeze(firing(:,l,k,1,1)))
%                         hold on
%                         plot(xaxis*180/pi,f_bl_mean,'r')
%                     end
%                 end
                [RFcenter,I] = max(firing(:));
                [i1,i2,RFcenterIdx,i3,i4] = ind2sub(size(firing),I);
                if plt == 1
                    axis([nana],[0 330 0 max(a(:,4))])
                    switch l
                        case 1
                            axes('Units','Normal');
                            h = title('Translation');
                            set(gca,'visible','off')
                            set(h,'visible','on')
%                             saveas(gca,[name,'_', num2str(ch),'_', num2str(unit) 't'],'png')
                        case 2
                            axes('Units','Normal');
                            h=title('Spirals');
                            set(gca,'visible','off')
                            set(h,'visible','on')
%                             saveas(gca,[name,'_', num2str(ch),'_', num2str(unit) 's'],'png')
                        case 3
                            axes('Units','Normal');
                            h= title('Shear');
                            set(gca,'visible','off')
                            set(h,'visible','on')
%                             saveas(gca,[name,'_', num2str(ch),'_', num2str(unit) 'd'],'png')
                    end
                end
            end
            for jjj = 1:rows
                responses1(:,:,jjj,:) = firing(:,:,(jjj-1)*columns+1:jjj*columns,j);
            end
            if ~(isempty(spikes)|| (strcmp(name(7),'c')&& microstimParadigm))% ~(isempty(spikes)|| isempty(spikes_bl)||strcmp(name(7),'c'))
                if ~RFmapping
                    figure
                    ShowSuperTuneAlt4(responses1,rows, columns,mean(firing_bl(:)),numMotion);
                    title([name, num2str(ch),' ', num2str(unit),', size: ',num2str(sizes(j))])
                    legend(['baseline = ',num2str(mean(firing_bl(:)))])
%                     saveas(gca,[name,'_', num2str(ch),'_', num2str(unit) 'TW'],'png')
                else
                    RF = squeeze(max(squeeze(responses1(:,1,:,:)),[],1));
                    x = -1*(columns-1)*sep/2:sep:(columns-1)*sep/2;
                    y = -1*(rows-1)*sep/2:sep:(rows-1)*sep/2;
                    figure
                    imagesc(RF)
                    title([name num2str(ch),' ', num2str(unit)])
                end
            end
        end
    end
    if sv == 1
        if microstim
            save(['E:\MT_MST\SuperTuneSpkTrains\test\',name(1:end-1),condition,num2str(ch),num2str(unit),'spktrain_bl.mat'],'spktrain_bl','Fs')
            save(['E:\MT_MST\SuperTuneSpkTrains\test\',name(1:end-1),condition,num2str(ch),num2str(unit),'spktrain.mat'],'spktrain','Fs','RFcenterIdx')
            save(['E:\MT_MST\SuperTuneFiringMatrix\test\',name(1:end-1),condition,num2str(ch),num2str(unit),'firingMat.mat'],'firing')
            save(['E:\MT_MST\SuperTuneFiringMatrix\test\',name(1:end-1),condition,num2str(ch),num2str(unit),'firingMat_bl.mat'],'firing_bl')
        else
            save(['E:\MT_MST\SuperTuneSpkTrains\ms_removed_all\',name,num2str(ch),num2str(unit),'spktrain_bl.mat'],'spktrain_bl','Fs')
            save(['E:\MT_MST\SuperTuneSpkTrains\ms_removed_all\',name,num2str(ch),num2str(unit),'spktrain.mat'],'spktrain','Fs','RFcenterIdx')
            save(['E:\MT_MST\SuperTuneFiringMatrix\ms_removed_all\',name,num2str(ch),num2str(unit),'firingMat.mat'],'firing')
            save(['E:\MT_MST\SuperTuneFiringMatrix\ms_removed_all\',name,num2str(ch),num2str(unit),'firingMat_bl.mat'],'firing_bl')
        end
    end
    clearvars -except plt sv v3categ flag nev DprimeV3 DprimeMT RFmapping sep sizeTuning firingAve channel firingGp chunum name Fs microstim condition rows columns microstimParadigm tstMS tendMS
end