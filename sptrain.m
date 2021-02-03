clear all; clc
close all
firingAve = zeros(32,2);
firingGp = zeros(32,9);

cd /mnt/F862AF7962AF3B70/MT_MST/MonkeyLab/HyperFlows/ytu332b/su/
f = dir('*.mat');
for i = 1:32
    fname(i).fn = f(i).name;
    fname(i).name = f(i).name;
    fname(i).name = erase(fname(i).name,'ytu332b');
    fname(i).name = erase(fname(i).name,'1N.mat');
end
index = [10;11;1;12;13;14;15;16;17;18;19;20;21;2;22;23;24;25;26;27;28;29;30;31;3;32;4;5;6;7;8;9];
for i = 1:32
    fn(index(i)).name = fname(i).fn; 
end
%%
for ch = 1
    name = 'ytu332b';
    unit = 1;
    NUMofMOTION = 3;
    load(fn(ch).name)
    featureIdx = 3; %the 8 directions for all motion types
    xaxis = unique(spikeMat(:,featureIdx));
    count = 1:length(xaxis);
    vec = zeros(length(unique(spikeMat(:,featureIdx))),1);
    gridIndeces = unique(spikeMat(:,7));
    numTrials = length(spikeMat(spikeMat(:,1)==-1000,1));
    startstim = 0;
    endstim = mean(stimLength)*1000;
    firing = zeros(length(xaxis),3,9);  
    firing_bl = firing;
    
    for l = 1:NUMofMOTION
        motionType = l;
        figure
        a = zeros(max(gridIndeces),4);
        for k = 1:max(gridIndeces)
            gridIndex = k;
            spikes = spikeMat(spikeMat(:,1)>startstim & spikeMat(:,4)==motionType & ...
                spikeMat(:,7)==gridIndex,:);
            trialsPerFeature=length(spikeMat(spikeMat(:,1)==-1000 & ...
                spikeMat(:,4)==motionType & spikeMat(:,7)==gridIndex,1))...
                /length(unique(spikeMat(:,3)));
          
            for i = count
                vec(i) = length(spikes(spikes(:,featureIdx)==xaxis(i),1));
                vec_rep_sep = spikes(spikes(:,featureIdx)==xaxis(i),:);
                spikes_bl = spikeMat(spikeMat(:,1)<=0 & spikeMat(:,1)>-1000 & spikeMat(:,4)==motionType &...
                    spikeMat(:,7)==gridIndex & spikeMat(:,featureIdx)==xaxis(i),:);
                firing_bl(i,l,k) = length(spikes_bl(:,1))/...
                    (trialsPerFeature*(mean(baseLineLength)));
                firing(i,l,k) = vec(i)/(trialsPerFeature*(endstim-startstim)/1000);
            end
            
            %%
            nana(k) = subplot(3,3,k);
            plot(nana(k),xaxis*180/pi,squeeze(firing(:,l,k)));
            hold on
            plot(xaxis*180/pi,firing_bl(:,l,k),'r')
            a(k,:) = axis();
            
        end
        axis([nana],[0 330 0 max(a(:,4))])
        switch l
            case 1
                axes('Units','Normal');
                h = title('Translation');
                set(gca,'visible','off')
                set(h,'visible','on')
            case 2
                axes('Units','Normal');
                h = title('Spirals');
                set(gca,'visible','off')
                set(h,'visible','on')
            case 3
                axes('Units','Normal');
                h = title('Shear');
                set(gca,'visible','off')
                set(h,'visible','on')
        end
    end
    %%
    metadata = 0; 
    responses1(:,:,1,:) = firing(:,:,1:3);
    responses1(:,:,2,:) = firing(:,:,4:6);
    responses1(:,:,3,:) = firing(:,:,7:9);
    firingAve(ch,:) = [mean(firing(:));mean(firing_bl(:))];
    firingGp(ch,:) = mean(squeeze(mean(firing,1)),1);
    if ~isempty(spikes)|| ~isempty(spikes_bl)
        figure
        ShowSuperTuneAlt4(responses1,metadata,mean(firing_bl(:)),NUMofMOTION);
        title([name num2str(ch),' ', num2str(unit)])
        legend(['baseline = ',num2str(mean(firing_bl(:)))])
    end
    
    m1(:,:) = firing(:,1,:); max_res1 = max(max(m1)); [loc{ch,1}(:,1) loc{ch,1}(:,2)] = find(m1 == max_res1); % Dir Pos
    m2(:,:) = firing(:,2,:); max_res2 = max(max(m2)); [loc{ch,2}(:,1) loc{ch,2}(:,2)] = find(m2 == max_res2);
    m3(:,:) = firing(:,3,:); max_res3 = max(max(m3)); [loc{ch,3}(:,1) loc{ch,3}(:,2)] = find(m3 == max_res3);
    clear m1 m2 m3 max_res
    
%     save(['D:\MT_MST\SuperTuneSpkTrains\mu',name,num2str(ch),num2str(unit),'spktrain_bl.mat'],'spktrain_bl','Fs')
%     save(['D:\MT_MST\SuperTuneSpkTrains\mu',name,num2str(ch),num2str(unit),'spktrain.mat'],'spktrain','Fs','RFcenterIdx')
%     save(['D:\MT_MST\SuperTuneFiringMatrix\mu',name,num2str(ch),num2str(unit),'firingMat.mat'],'firing')
%     save(['/mnt/F862AF7962AF3B70/MT_MST/STFM/',name,num2str(ch),num2str(unit),'firingMat.mat'],'firing')

%     clearvars -except firingAve channel firingGp
end

