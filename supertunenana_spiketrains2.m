clear all; close all; clc
firingAve=zeros(32,2);
firingGp=zeros(32,9);
Fs=10000;
microstim=0;
condition='S';
name='ytu310a';
load(['E:\MT_MST\Synchrony things\chunum\',name(1:end-1),'Chunum.mat'])
% ci2=1;
chunum = chunum;%unit3 3 5 6 14 15 18 19 20 26 27 29 unit4 18 20 27
%chunum=chunum(goodcells,:);
load(['E:\MT_MST\Plexon\RFiles\',name,'_TrialStructure.mat'])
rows = file.taskDialogValues.superTuneRows;
columns = file.taskDialogValues.superTuneColumns;

for ci=1:length(chunum)
    counting=1;
    ch=chunum(ci,1);
    unit=chunum(ci,2);
    if microstim
        load([name(1:end-1),condition,num2str(ch),num2str(unit),'N.mat']);
        if strcmp(condition,'N')
            spikeMat=spikeMatnostim;
        else
            spikeMat=spikeMatstim;
        end
        stimInfo=load(['E:\MT_MST\Plexon\RFiles\',name,num2str(ch),num2str(unit),'N.mat']);
        stimLength=stimInfo.stimLength;
        baseLineLength=stimInfo.baseLineLength;
    else
        [spikeMat,numspikes,stimLength,baseLineLength]=plx2NFile_superTune(name,...
            ch,unit,1);
    end
    % spikeMat2(:,2)=spikeMat2(:,2)+spikeMat1(end,2);
    % spikeMat=[spikeMat1;spikeMat2];
    % numspikes=numspikes1+numspikes2;
    % stimLength=[stimLength1,stimLength2];
    % baseLineLength=[baseLineLength1,baseLineLength2];
    t1=0;
    t2=mean(stimLength);
    time=floor((t2-t1)*Fs);
    time_bl=floor(mean(baseLineLength)*Fs);
    %in spikemat column 2 trial number, column 3 direction of motion (8
    %directions), column 4 type of motion (3 types), column 5 speed of
    %translation,
    %spike times are aligned from stim on. 0 means stim on
    featureIdx=3; %the 8 directions for all motion types
    xaxis=unique(spikeMat(:,featureIdx));
    numMotion=max(unique(spikeMat(:,4)));
    count=1:length(xaxis);
    vec=zeros(length(unique(spikeMat(:,featureIdx))),1);
    gridIndeces=unique(spikeMat(:,7));
    numTrials=length(spikeMat(spikeMat(:,1)==-1000,1));
    startstim=0;
    endstim=mean(stimLength)*1000;
    %motionType=2;
    %gridIndex=5;
    %trialsPerFeature=length(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==1
    %& spikeMat(:,7)==1,1))/length(unique(spikeMat(:,3)));
    
    firing=zeros(length(xaxis),numMotion,rows*columns);
    firing_bl=firing;
    %firing2=zeros(8,3,9,trialsPerFeature);
    % spikes_bl=spikeMat(spikeMat(:,1)<=0 ,:);
    % firing_bl=length(spikes_bl(:,1))/(numTrials*(mean(baseLineLength)));
    for l=1:numMotion
        motionType=l;
        figure
        a=zeros(max(gridIndeces),4);
        for k=1:max(gridIndeces)
            gridIndex=k;
            
            spikes=spikeMat(spikeMat(:,1)>startstim & spikeMat(:,4)==motionType & ...
                spikeMat(:,7)==gridIndex,:);
            trialsPerFeature=floor(length(unique(spikeMat(spikeMat(:,1)==-1000 & ...
                spikeMat(:,4)==motionType & spikeMat(:,7)==gridIndex,2)))...
                /length(unique(spikeMat(:,3))));
            if counting==1
                spktrain=zeros(time,length(xaxis),numMotion,rows*columns,trialsPerFeature);
                spktrain_bl=zeros(time_bl,length(xaxis),numMotion,rows*columns,trialsPerFeature);
                counting=2;
            end
            for i=count
                % hist(spikeMat(spikeMat(:,3)==i,1));
                vec(i)=length(spikes(spikes(:,featureIdx)==xaxis(i),1));
                vec_rep_sep=spikes(spikes(:,featureIdx)==xaxis(i),:);
                spikes_bl=spikeMat(spikeMat(:,1)<=0 & spikeMat(:,1)>=(0-time_bl)*1000/Fs & spikeMat(:,4)==motionType &...
                    spikeMat(:,7)==gridIndex & spikeMat(:,featureIdx)==xaxis(i),:);
                trialsPerThisFeature=floor(length(unique(spikeMat(spikeMat(:,1)==-1000 & ...
                    spikeMat(:,4)==motionType & spikeMat(:,7)==gridIndex  & spikeMat(:,featureIdx)==xaxis(i),2))));
                firing_bl(i,l,k)=length(spikes_bl(:,1))/...
                    (trialsPerThisFeature*(mean(baseLineLength)));
                for o=1:trialsPerThisFeature
                    trialIdx=unique(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==motionType &...
                        spikeMat(:,7)==gridIndex & spikeMat(:,featureIdx)==xaxis(i),2));
                    spktrain1=zeros(time,1);
                    indices=round(vec_rep_sep(vec_rep_sep(:,2)==trialIdx(o),1)*(Fs/1000));
                    indices=indices(indices<=time);
                    spktrain1(indices)=1;
                    spktrain(:,i,l,k,o)=spktrain1;
                    spktrain2=zeros(time_bl,1);
                    if trialIdx(o)> length(baseLineLength)
                        indices_bl=round((spikes_bl(spikes_bl(:,2)==trialIdx(o),1)+(mean(baseLineLength)*1000))*(Fs/1000));
                    else
                        indices_bl=round((spikes_bl(spikes_bl(:,2)==trialIdx(o),1)+(baseLineLength(trialIdx(o))*1000))*(Fs/1000));
                    end
                    indices_bl=indices_bl((indices_bl>0)&(indices_bl<=time_bl));
                    spktrain2(indices_bl)=1;
                    spktrain_bl(:,i,l,k,o)=spktrain2;
                end
                firing(i,l,k)=vec(i)/(trialsPerThisFeature*(endstim-startstim)/1000);
                %firing_bl(i)=length(spikes_bl(spikes_bl(:,featureIdx)==xaxis(i),1))/(trialsPerFeature*startstim/1000);
                %ylim([0 100])
            end
            
            
            %%
            nana(k)=subplot(rows,columns,k); %was 3 l k
            plot(nana(k),xaxis*180/pi,squeeze(firing(:,l,k)));
            hold on
            plot(xaxis*180/pi,firing_bl(:,l,k),'r')
            %polar(xaxis,vec)
            a(k,:)=axis();
            
        end
        [RFcenter,I]=max(firing(:));
        [i1,i2,RFcenterIdx]= ind2sub(size(firing),I);
        axis([nana],[0 330 0 max(a(:,4))])
        switch l
            case 1
                axes('Units','Normal');
                h = title('Translation');
                set(gca,'visible','off')
                set(h,'visible','on')
            case 2
                axes('Units','Normal');
                h=title('Spirals');
                set(gca,'visible','off')
                set(h,'visible','on')
            case 3
                axes('Units','Normal');
                h= title('Shear');
                set(gca,'visible','off')
                set(h,'visible','on')
        end
    end
    % save(['C:\research\data\SuperTuneSpkTrains\',name,num2str(chunum(ci,1)),num2str(chunum(ci,2)),'spktrain_bl.mat'],'spktrain_bl')
    % save(['C:\research\data\SuperTuneSpkTrains\',name,num2str(chunum(ci,1)),num2str(chunum(ci,2)),'spktrain.mat'],'spktrain')
    if microstim
        save(['E:\MT_MST\SuperTuneSpkTrains\',name(1:end-1),condition,num2str(ch),num2str(unit),'spktrain_bl.mat'],'spktrain_bl','Fs')
        save(['E:\MT_MST\SuperTuneSpkTrains\',name(1:end-1),condition,num2str(ch),num2str(unit),'spktrain.mat'],'spktrain','Fs','RFcenterIdx')
        save(['E:\MT_MST\SuperTuneFiringMatrix\',name(1:end-1),condition,num2str(ch),num2str(unit),'firingMat.mat'],'firing')
    else
        save(['E:\MT_MST\SuperTuneSpkTrains\',name,num2str(ch),num2str(unit),'spktrain_bl.mat'],'spktrain_bl','Fs')
        save(['E:\MT_MST\SuperTuneSpkTrains\',name,num2str(ch),num2str(unit),'spktrain.mat'],'spktrain','Fs','RFcenterIdx')
        save(['E:\MT_MST\SuperTuneFiringMatrix\',name,num2str(ch),num2str(unit),'firingMat.mat'],'firing')
    end
    %%
    metadata=0;
    %firing_bl=0.1;
    %responses=reshape(firing,[8 3 3 3]);
    % if strcmp(name,'ytu196a')
    % responses1(:,:,1,:)=firing(1:4,:,1:3);
    % responses1(:,:,2,:)=firing(1:4,:,4:6);
    % responses1(:,:,3,:)=firing(1:4,:,7:9);
    % else
    
    responses1(:,:,1,:)=firing(:,:,1:3);
    responses1(:,:,2,:)=firing(:,:,4:6);
    responses1(:,:,3,:)=firing(:,:,7:9);
    firingAve(ch,:)=[mean(firing(:));mean(firing_bl(:))];
    firingGp(ch,:)=mean(squeeze(mean(firing,1)),1);
    %end
    if ~isempty(spikes)|| ~isempty(spikes_bl)
        figure
        ShowSuperTuneAlt4(responses1,rows,columns,mean(firing_bl(:)),numMotion);
        title([name num2str(ch),' ', num2str(unit)])
        legend(['baseline = ',num2str(mean(firing_bl(:)))])
    end
    clearvars -except firingAve channel firingGp chunum name Fs microstim condition rows columns
end
mana=[firingAve,zeros(32,1),firingGp];
figure
scatter(firingAve(:,1),firingAve(:,2))
xlabel('stim on')
ylabel('baseline')
close all

%clear variables
%end
%savefig([name num2str(ch) num2str(unit) '_TuningWheels.fig'])
%end
%end
%axis([nana],[min(min(a(:,1))) max(max(a(:,2))) min(min(a(:,3))) max(max(a(:,4)))])
%axis([min(a(:,1)) max(a(:,2)) min(a(:,3)) max(a(:,4))])