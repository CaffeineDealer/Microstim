%% Z-score normalization
close all
clear all; clc
Fs = 10000;
name = 'ytu337';
path = 'E:\MT_MST\SuperTuneFiringMatrix\ms_rmvd_150300\';
pathspk = 'E:\MT_MST\SuperTuneSpkTrains\ms_rmvd_150300\';
% clib = 'E:\MT_MST\Microstim\Cell_Lib\Candidate Cells\';
% ch = struct2array(load([clib,'CLib_',name,'.mat']));
FsMS = 400; 
nPulse = 5;
PulseDur = (1/FsMS) * nPulse * Fs; 
OnsetMS = 0.150;
tstMS = OnsetMS * Fs;
tendMS = tstMS + PulseDur;
send = 0.3 * Fs;
ch = 1:64;
% load('E:\MT_MST\Plexon\RFiles\nevlib.mat')

for i = 1:length(ch)
    a = struct2array(load([path,name,'a',num2str(ch(i)),'1','firingMat.mat'])); a = a(:,1:2,:);
    b = struct2array(load([path,name,'b',num2str(ch(i)),'1','firingMat.mat'])); b = b(:,1:2,:);
    
    a_bl = struct2array(load([path,name,'a',num2str(ch(i)),'1','firingMat_bl.mat'])); a_bl = a_bl(:,1:2,:);
    b_bl = struct2array(load([path,name,'b',num2str(ch(i)),'1','firingMat_bl.mat'])); b_bl = b_bl(:,1:2,:);
    
    load([pathspk,name,'a',num2str(ch(i)),'1','spktrain.mat']); aspk = spktrain(tendMS:send,:,1:2,:,:);
    load([pathspk,name,'b',num2str(ch(i)),'1','spktrain.mat']); bspk = spktrain(tendMS:send,:,1:2,:,:);
    
    load([pathspk,name,'a',num2str(ch(i)),'1','spktrain_bl.mat']); aspk_bl = spktrain_bl(:,:,1:2,:,:);
    load([pathspk,name,'b',num2str(ch(i)),'1','spktrain_bl.mat']); bspk_bl = spktrain_bl(:,:,1:2,:,:);
    
    a = a - mean(a_bl(:));
    a_bl = a_bl - mean(a_bl(:));
    za = a; za(za<0) = 0;
    za_bl = a_bl; za_bl(za_bl<0) = 0;
    
    b = b - mean(b_bl(:));
    b_bl = b_bl - mean(b_bl(:));
    zb = b; zb(zb<0) = 0;
    zb_bl = b_bl; zb_bl(zb_bl<0) = 0;
    
%     za = z_score(a,2); za_bl = z_score(a_bl,2);
%     zb = z_score(b,2); zb_bl = z_score(b_bl,2);

    [NMS MS] = MS_NMS(za,zb,8,2,9);
    [NMS_bl MS_bl] = MS_NMS(za_bl,zb_bl,8,2,9);
    
    [NMSspk MSspk] = loadspk(aspk,bspk,8,2,9,'stim');
    [NMSspk_bl MSspk_bl] = loadspk(aspk_bl,bspk_bl,8,2,9,'bl');
    
    firing = NMS;
    firing_bl = NMS_bl;
    spktrain = NMSspk;
    spktrain_bl = NMSspk_bl; 
    save([path,name,'N',num2str(ch(i)),num2str(1),'firingMat.mat'],'firing')
    save([path,name,'N',num2str(ch(i)),num2str(1),'firingMat_bl.mat'],'firing_bl')
    save([pathspk,name,'N',num2str(ch(i)),num2str(1),'spktrain.mat'],'spktrain')
    save([pathspk,name,'N',num2str(ch(i)),num2str(1),'spktrain_bl.mat'],'spktrain_bl')
    % ttest to keep/discard cells
    baseline = squeeze(sum(spktrain_bl,1))*Fs/size(spktrain_bl,1);
    allstimfir = squeeze(sum(spktrain,1))*Fs/size(spktrain,1);
    [h,p] = ttest(baseline(:),allstimfir(:));
    keepCriteria(i,1) = (p <= 0.05);
    
    clear firing firing_bl spktrain spktrain_bl
    
    firing = MS;
    firing_bl = MS_bl;
    spktrain = MSspk;
    spktrain_bl = MSspk_bl;
    save([path,name,'S',num2str(ch(i)),'1','firingMat.mat'],'firing')
    save([path,name,'S',num2str(ch(i)),num2str(1),'firingMat_bl.mat'],'firing_bl')
    save([pathspk,name,'S',num2str(ch(i)),num2str(1),'spktrain.mat'],'spktrain')
    save([pathspk,name,'S',num2str(ch(i)),num2str(1),'spktrain_bl.mat'],'spktrain_bl')
    % ttest to keep/discard cells
    baseline = squeeze(sum(spktrain_bl,1))*Fs/size(spktrain_bl,1);
    allstimfir = squeeze(sum(spktrain,1))*Fs/size(spktrain,1);
%     [h,p] = ttest(baseline(:),allstimfir(:));
%     keepCriteria(i,1) = (p <= 0.05);
    
    bl = [NMS_bl MS_bl];
    mbl = mean(bl(:));
    tcplot([0:45:315],NMS,1,3,3,'Normal',mbl,ch(i),keepCriteria(i,1))%ch(i) - nevlib(i,2)
    tcplot([0:45:315],NMS,2,3,3,'Normal',mbl,ch(i),keepCriteria(i,1))
end
