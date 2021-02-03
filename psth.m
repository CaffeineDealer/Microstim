clear all; clc
path = 'E:\MT_MST\SuperTuneSpkTrains\';
path2 = 'E:\MT_MST\SuperTuneFiringMatrix\';
name = 'ytu309a';
condition = 'S';
load(['E:\MT_MST\Synchrony things\chunum\',name(1:end-1),'Chunum.mat'])
bin = 100; %10 ms
stimCh = 52;
u = 1;
stimfiring = load([path2,name,num2str(stimCh),num2str(u),'firingMat']);
firingst = stimfiring.firing;
%%
[maxfir,I]=max(firingst(:));
[pdir,pmot,RFcenterIdx]= ind2sub(size(firingst),I);
numdir=size(firingst,1);
if pdir>numdir/2
    npdir=pdir-(numdir/2);
else
    npdir=pdir+(numdir/2);
end
%%
figure
count2=1;
for range = 1:length(chunum(:,1))    
    count=1;
    count3=1;
    ch=chunum(range,1);
    unit=1;%chunum(range,2);
    load([path,name(1:end-1),'a',num2str(ch),num2str(unit),'spktrain.mat']);
    soso=reshape(spktrain,[size(spktrain,1) size(spktrain,2)*size(spktrain,3)*...
        size(spktrain,4)*size(spktrain,5)]);
    sosoP=reshape(spktrain(:,pdir,pmot,:,:),[size(spktrain,1)  size(spktrain,4)*size(spktrain,5)]);
    sosoNP=reshape(spktrain(:,npdir,pmot,:,:),[size(spktrain,1)  size(spktrain,4)*size(spktrain,5)]);
    load([path,name(1:end-1),'N',num2str(ch),num2str(unit),'spktrain.mat']);
    nono=reshape(spktrain,[size(spktrain,1) size(spktrain,2)*size(spktrain,3)*...
        size(spktrain,4)*size(spktrain,5)]);
    nonoP=reshape(spktrain(:,pdir,pmot,:,:),[size(spktrain,1)  size(spktrain,4)*size(spktrain,5)]);
    nonoNP=reshape(spktrain(:,npdir,pmot,:,:),[size(spktrain,1)  size(spktrain,4)*size(spktrain,5)]);   
    for k=1:bin:size(soso,1)-bin
        lala=sum(soso((k:k+bin),:),1);
        psthst(count2,count) = mean(lala);
        lalaP=sum(sosoP((k:k+bin),:),1);
        psthstP(count2,count)=mean(lalaP);
        lalaNP=sum(sosoNP((k:k+bin),:),1);
        psthstNP(count2,count)=mean(lalaNP);
        tata=sum(nono((k:k+bin),:),1);
        psthn(count2,count)=mean(tata);
        tataP=sum(nonoP((k:k+bin),:),1);
        psthnP(count2,count)=mean(tataP);
        tataNP=sum(nonoNP((k:k+bin),:),1);
        psthnNP(count2,count)=mean(tataNP);
        count=count+1;
    end
    plot(psthstP(count2,:)); hold on
    title([name(1:end-1), num2str(ch),' ', num2str(unit)])
end
