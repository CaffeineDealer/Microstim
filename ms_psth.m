clear all; clc
path = 'E:\MT_MST\SuperTuneSpkTrains\';
path2 = 'E:\MT_MST\SuperTuneFiringMatrix\';
name = 'ytu310S';
% stimCh = 31;
% 2,52,54,33,56,57,57,52,58,50,60,43,51,64,16
u = 1;
mtype = 1; 
% load(['E:\MT_MST\Synchrony things\chunum\',name(1:end-1),'Chunum.mat'])
load(['E:\MT_MST\Microstim\PSTH\',name(1:end-1),'Chunum.mat'])
% stimfiring = load([path2,name(1:end-1),'a',num2str(stimCh),num2str(u),'firingMat']);
% firingst = stimfiring.firing;
bin = 100; 
% [pdir,pmot,ppos,npdir] = preferred(firingst);
%%
chu = chunum(chunum(:,2) == 1);
nch = max(chunum(:));
for range = 1:length(chu)    
    ch = chu(range,1);
    load([path,name(1:end-1),'S',num2str(ch),num2str(u),'spktrain.mat']);
    load([path,name(1:end-1),'S',num2str(ch),num2str(u),'spktrain_bl.mat'])
    stimP = reshape(spktrain(:,:,:,:,:),[size(spktrain,1)  size(spktrain,2)*size(spktrain,3)*size(spktrain,4)*size(spktrain,5)]);
    stimP_bl = reshape(spktrain_bl(:,:,:,:,:),[size(spktrain_bl,1)  size(spktrain_bl,2)*size(spktrain_bl,3)*size(spktrain_bl,4)*size(spktrain_bl,5)]);
%     stimNP = reshape(spktrain(:,:,:,:,:),[size(spktrain,1)  size(spktrain,2)*size(spktrain,3)*size(spktrain,4)*size(spktrain,5)]);
%     stimP = reshape(spktrain(:,pdir,pmot,:,:),[size(spktrain,1)  size(spktrain,4)*size(spktrain,5)]);
%     stimNP = reshape(spktrain(:,npdir,pmot,:,:),[size(spktrain,1)  size(spktrain,4)*size(spktrain,5)]);    
    clear spktrain spktrain_bl
    load([path,name(1:end-1),'N',num2str(ch),num2str(u),'spktrain.mat']);
    load([path,name(1:end-1),'N',num2str(ch),num2str(u),'spktrain_bl.mat']);
    nstimP = reshape(spktrain(:,:,:,:,:),[size(spktrain,1)  size(spktrain,2)*size(spktrain,3)*size(spktrain,4)*size(spktrain,5)]);
    nstimP_bl = reshape(spktrain_bl(:,:,:,:,:),[size(spktrain_bl,1)  size(spktrain_bl,2)*size(spktrain_bl,3)*size(spktrain_bl,4)*size(spktrain_bl,5)]);
%     nstimNP = reshape(spktrain(:,:,:,:,:),[size(spktrain,1)  size(spktrain,2)*size(spktrain,3)*size(spktrain,4)*size(spktrain,5)]);
%     nstimP = reshape(spktrain(:,pdir,pmot,:,:),[size(spktrain,1)  size(spktrain,4)*size(spktrain,5)]);
%     nstimNP = reshape(spktrain(:,npdir,pmot,:,:),[size(spktrain,1)  size(spktrain,4)*size(spktrain,5)]);
    for i = 1:size(stimP,2)
        SP(ch,:,i) = psth_sp(stimP,bin,i);
        NSP(ch,:,i) = psth_sp(nstimP,bin,i);
        SP_bl(ch,:,i) = psth_sp(stimP_bl,bin,i);
        NSP_bl(ch,:,i) = psth_sp(nstimP_bl,bin,i);
%         SNP(ch,:,i) = psth_sp(stimNP,bin,i);
%         NSNP(ch,:,i) = psth_sp(nstimNP,bin,i);
    end
end
%%
sp_m = zeros(2,size(SP,1),size(SP_bl,2)+size(SP,2));
for i = 1:size(SP,1)
    sp_m(1,i,1:size(SP_bl,2)) = (mean(squeeze(SP_bl(i,:,:)),2));
    sp_m(1,i,size(SP_bl,2)+1:end) = (mean(squeeze(SP(i,:,:)),2));
    sp_m(2,i,1:size(NSP_bl,2)) = (mean(squeeze(NSP_bl(i,:,:)),2));
    sp_m(2,i,size(NSP_bl,2)+1:end) = (mean(squeeze(NSP(i,:,:)),2));
%     sp_m(3,i,:) = (mean(squeeze(NSP(i,:,:)),2));
%     sp_m(4,i,:) = (mean(squeeze(NSNP(i,:,:)),2));
end
%%
% time = 0:bin/Fs:size(spktrain,1)/Fs;
% figure
% for i = 1:size(SP,1)
%     ctl = (mean(squeeze(SP(i,:,:)),2));
%     plot(time,ctl,'LineWidth',1.5); hold on
% end
% xlabel 'Time(ms)'
% ylabel 'Mean Spike Count'
%%
figure
time1 = -size(spktrain_bl,1)/Fs:bin/Fs:0;
time2 = 0:bin/Fs:size(spktrain,1)/Fs;
time = [time1 time2];
STime = find(time == 0.15);
Ston = find(time == 0);
dict = {'Stim - Microstim','Stim - No Microstim'};
for i = 1:2
    subplot(2,1,i)
    imagesc(squeeze(sp_m(i,:,:))); colorbar; hold on
    line([Ston Ston],[0 max(chu)],'LineWidth',2,'Color','g')
    line([STime STime],[0 max(chu)],'LineWidth',2,'Color','r')
    set(gca,'YDir','normal')
    title(sprintf('%s - %s',name(1:end-1),dict{i}))
    xlabel 'Time'
    ylabel 'Ch'
    colorbar
end
