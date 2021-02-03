function [stim stim_bl sz_spk sz_spk_bl] = spkt(path,name,condition,ch,u,mtype,pos)


load([path,name(1:end-1),condition,num2str(ch),num2str(u),'spktrain.mat'])
load([path,name(1:end-1),condition,num2str(ch),num2str(u),'spktrain_bl.mat'])
stim = reshape(spktrain(:,:,mtype,pos,:),[size(spktrain,1)  size(spktrain,2)*size(spktrain,5)]);
stim_bl = reshape(spktrain_bl(:,:,mtype,pos,:),[size(spktrain_bl,1)  size(spktrain_bl,2)*size(spktrain_bl,5)]);

sz_spk = size(spktrain,1);
sz_spk_bl = size(spktrain_bl,1);