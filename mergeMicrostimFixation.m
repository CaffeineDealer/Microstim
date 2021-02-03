clear all; clc
file = 'ytu308c';
unit = 1;
load(['E:\MT_MST\Synchrony things\chunum\',file(1:end-1),'Chunum.mat'])
%%
for mm = 1:length(chunum)
    ch = chunum(mm,1);
    unit = 1;%chunum(mm,2);
    spikeMatt = load(['E:\MT_MST\MonkeyLab\HyperFlows\',file,'\',file,num2str(ch),num2str(unit),'N.mat']);
    spikeMat1 = spikeMatt.spikeMat;
    
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
        spikeMatnostim1 = spikeMat1(spikeMat1(:,3) == nonstimidx(i),:);
        spikeMatnostim = [spikeMatnostim;spikeMatnostim1];
        spikeMatstim1 = spikeMat1(spikeMat1(:,3) == stimidx(i),:);
        spikeMatstim = [spikeMatstim;spikeMatstim1];
    end
    spikeMatnostim = sortrows(spikeMatnostim,2);
    spikeMatstim = sortrows(spikeMatstim,2);
    trialnum = 0;
    trialnum2 = 0;
    for k = 1:length(spikeMatnostim(:,1))
        if spikeMatnostim(spikeMatnostim(k,1) == -1000)
            trialnum = trialnum + 1;
            spikeMatnostim(k,2) = trialnum;
        else
            spikeMatnostim(k,2) = trialnum;
        end
    end
    for l = 1:length(spikeMatstim(:,1))
        if spikeMatstim(spikeMatstim(l,1) ==- 1000)
            trialnum2 = trialnum2 + 1;
            spikeMatstim(l,2) = trialnum2;
        else
            spikeMatstim(l,2) = trialnum2;
        end
        
    end
    save(['E:\MT_MST\Plexon\RFiles\',file(1:end-1),'Sf',num2str(ch),num2str(unit),'N.mat'],'spikeMatstim');
    save(['E:\MT_MST\Plexon\RFiles\',file(1:end-1),'Nf',num2str(ch),num2str(unit),'N.mat'],'spikeMatnostim');
end
