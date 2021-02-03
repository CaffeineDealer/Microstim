clear all; clc
% takes 2 interleaved supertune spikeMat and spike trains and creates non
% stim spiketrain and spike mat and another only stim
%feature for microstim in file 1, feature for microstim in file 2
% for spike mat file 1 odd numbered direction indeces
file1='ytu310a';
file2='ytu310b';
% chunum = [1:32]';
% switch file1 
%     case 'ytu310a'
%         chunum = [51:58]';
%     case 'ytu316a'
%         chunum = [52;54;56;57;58;59;60];
%     case 'ytu321a'
%         chunum = [54:60]';
%     case 'ytu323a'
%         chunum = [49:55]';
%     case 'ytu333a'
%         chunum = [48:55]';
%     case 'ytu335a'
%         chunum = [59:64]';
% end
load(['E:\MT_MST\Microstim\Cell_Lib\Candidate Cells\','CLib_',file1(1:end-1),'.mat'])
chunum = CLib;
unit=1;
% load(['E:\MT_MST\Synchrony things\chunum\',file1(1:end-1),'Chunum.mat'])
% chunum = [59:64]';
% load(['/home/caffeinedealer/Documents/MicroStim/Synchrony things/chunum/',file1(1:end-1),'Chunum.mat'])
for mm = 1:length(chunum)
    ch = chunum(mm,1);
    unit = 1;%chunum(mm,2);
    spikeMatt = load(['E:\MT_MST\MonkeyLab\HyperFlows\',file1,'\',file1,num2str(ch),num2str(unit),'N.mat']);
    spikeMat1 = spikeMatt.spikeMat;
    spikeMats = load(['E:\MT_MST\MonkeyLab\HyperFlows\',file2,'\',file2,num2str(ch),num2str(unit),'N.mat']);
    spikeMat2 = spikeMats.spikeMat;
    angles = sort(unique(spikeMat1(:,3)));
    numangles = length(angles);
    if numangles == 4
        stimidx = angles(logical([1 0 1 0]));%angles(1,3);
        nonstimidx = angles(~logical([1 0 1 0]));%angles(2,4);
    elseif numangles == 6
        stimidx = angles(logical([1 0 1 0 1 0]));%angles(1,3,5);
        nonstimidx = angles(~logical([1 0 1 0 1 0]));%angles(2,4,6);
    elseif numangles == 8
        stimidx = angles(logical([1 0 1 0 1 0 1 0]));%angles(1,3,5,7);
        nonstimidx = angles(~logical([1 0 1 0 1 0 1 0]));%angles(2,4,6,8);
    end
    spikeMatnostim = [];
    spikeMatstim = [];
    for i = 1:length(stimidx)
        spikeMatnostim1 = spikeMat1(spikeMat1(:,3)==nonstimidx(i),:);
        spikeMatnostim2 = spikeMat2(spikeMat2(:,3)==stimidx(i),:);
        spikeMatnostim = [spikeMatnostim;spikeMatnostim1;spikeMatnostim2];
        spikeMatstim1 = spikeMat1(spikeMat1(:,3)==stimidx(i),:);
        spikeMatstim2 = spikeMat2(spikeMat2(:,3)==nonstimidx(i),:);
        spikeMatstim = [spikeMatstim;spikeMatstim1;spikeMatstim2];
    end
    spikeMatnostim = sortrows(spikeMatnostim,2);
    spikeMatstim = sortrows(spikeMatstim,2);
    trialnum = 0;%1:length(spikeMat1(spikeMat1(:,1)==-1000));
    trialnum2 = 0;
    for k = 1:length(spikeMatnostim(:,1))
        if spikeMatnostim(spikeMatnostim(k,1)==-1000)
            trialnum = trialnum+1;
            spikeMatnostim(k,2) = trialnum;
        else
            spikeMatnostim(k,2) = trialnum;
        end
    end
    for l = 1:length(spikeMatstim(:,1))
        if spikeMatstim(spikeMatstim(l,1)==-1000)
            trialnum2=trialnum2+1;
            spikeMatstim(l,2)=trialnum2;
        else
            spikeMatstim(l,2)=trialnum2;
        end
        
    end
    save(['E:\MT_MST\Plexon\RFiles\',file1(1:end-1),'S',num2str(ch),num2str(unit),'N.mat'],'spikeMatstim');
    save(['E:\MT_MST\Plexon\RFiles\',file2(1:end-1),'N',num2str(ch),num2str(unit),'N.mat'],'spikeMatnostim');
end
