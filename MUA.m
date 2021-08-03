clear all; clear cache; close all; clc
cd E:\MT_MST\Plexon\PLXfiles
%% Merge a & b + create spikewaveforms
name = 'ytu310';
a = load(sprintf('%sa.mat',name));
b = load(sprintf('%sb.mat',name));
Fs = a.params.Fsds;

allcomb = a.params.comb.allcomb;
if size(allcomb,1) == 144 
    NMSa = [19,36;55,72;91,108;127,144];
    NMSb = [1,18;37,54;73,90;109,126];
    counter = 18;
elseif size(allcomb,1) == 216
    NMSa = [28,54;82,108;136,162;190,216];
    NMSb = [1,27;55,81;109,135;163,189];
    counter = 27;
end
Xa = {a.xsorted.xfs}';
Xb = {b.xsorted.xfs}';

spkt = [];
for i = 1:4
    r1 = NMSb(i,1);
    r2 = NMSb(i,2);
    r3 = NMSa(i,1);
    for j = r1:r2
        spkt(j).NMS = Xb{r1,1};
        spkt(j+counter).NMS = Xa{r3,1};
        spkt(j).MS = Xa{r1,1};
        spkt(j+counter).MS = Xb{r3,1};
        r1 = r1 + 1;
        r3 = r3 + 1;
    end
end
%% Covering missing trials with zeros
for i = 1:size(spkt,2)
    if size(spkt(i).NMS,3) ~= 4
        missC = 4 - size(spkt(i).NMS,3);
        spkt(i).NMS(:,:,end + missC) = zeros;
    end
    if size(spkt(i).MS,3) ~= 4
        missC = 4 - size(spkt(i).MS,3);
        spkt(i).MS(:,:,end + missC) = zeros;
    end
end
%% Creating a 6D matrix ch * npoints * dir * motion * position * repetition
NMS = [];
MS = [];
NMScell = {spkt.NMS}';
MScell = {spkt.MS}';
i = 1;
for j = 1:8 % #direction 
    for z = 1:size(unique(allcomb(:,2)),1) % #motion
        for w = 1:size(unique(allcomb(:,3)),1) % #position
            NMS(j,z,w,:,:,:) = double(NMScell{i,1});
            MS(j,z,w,:,:,:) = double(MScell{i,1});
            i = i + 1;
        end
    end
end
NMS = permute(NMS,[4 5 1 2 3 6]);
MS = permute(MS,[4 5 1 2 3 6]);
NMS(:,4001:5000,:,:,:,:) = [];
MS(:,4001:5000,:,:,:,:) = [];
%% remove microstim artifact
lbthr = 0.0005;
zeroSize = 0.03 * Fs; % 30 msec 
blackoutNMS = findblackout(NMS,lbthr,zeroSize);
blackoutMS = findblackout(MS,lbthr,zeroSize);
%% Threshold for spike sorting
stdmin = 3; % min threshold for detection
stdmax = 10; % max threshold for detection
[thrnms thrmaxnms] = findthreshold(NMS,stdmin,stdmax);
[thrms thrmaxms] = findthreshold(MS,stdmin,stdmax);
%% Generate spike trains
spktNMS = [];
spktMS = [];
nrept = 4;
spktNMS = findspkt(NMS,nrept,thrnms,thrmaxnms);%thrnms
spktMS = findspkt(MS,nrept,thrms,thrmaxms);
%% Firing Rate and Discarding blackouts
[frbl frst prCorrect fr_bl] = FireRate_bin(spktNMS,blackoutNMS,Fs);
save(['E:\MT_MST\Microstim\MST-MUA\',sprintf('%sN.mat',name)],'frst','frbl','spktNMS','blackoutNMS','prCorrect')
frbl = []; frst = []; prCorrect = []; fr_bl = [];
[frbl frst prCorrect fr_bl] = FireRate_bin(spktMS,blackoutMS,Fs);
save(['E:\MT_MST\Microstim\MST-MUA\',sprintf('%sS.mat',name)],'frst','frbl','spktNMS','blackoutNMS','prCorrect')
soundsc(rand(3000,1))
return
% [frbl frst prCorrect] = FireRate(spktNMS,blackoutNMS,Fs);
%% Plot Tuning Curves
clc
dir = [0:45:315];
for i = 1:32
    tcplotMUA(dir,frst,1,3,3,frbl,i,prCorrect)
    tcplotMUA(dir,frst,2,3,3,frbl,i,prCorrect)
    % T-Test to find good cell
%     pval(i,1) = GoodCellTtest(frst,fr_bl,i);
end
%% Generate Random firing rates to test algorithm 
frst = [];
frbl = [];
x = zeros(64,8000,8,3,9,4);
for i = 1:64
    ct = floor(rand(1)*20 + 20);
    a = 4500;
    b = 5500;
    r = floor((b-a).*rand(ct,1) + a);
    x(i,r,2,2,2,1) = 1;
end
[frbl frst prCorrect] = FireRate_bin(x,bout,Fs);
% [frbl frst prCorrect] = FireRate(x,bout,Fs);
%%
% spktNMS = x;
% blackoutNMS = bout;
ch = 51;
mo = 2;
pos = 1;
tr = 1;
figure
for i = 1:8 
    subplot(4,2,i)
    if blackoutNMS(ch,i,mo,pos,tr) == 1
        ctl = 'fail';
        plot((NMS(ch,:,i,mo,pos,tr)),'r')
%         plot(spktNMS(ch,:,i,mo,pos,tr),'r')
    elseif blackoutNMS(ch,i,mo,pos,tr) == 0
        ctl = 'pass';
%         plot((NMS(ch,:,i,mo,pos,tr)),'b')
        plot(spktNMS(ch,:,i,mo,pos,tr),'b')
    end
    title(sprintf('%s',ctl))
end
%%
dir = 3;
mo = 2;
pos = 1;
tr = 1;
cc = [1:64];
figure
for i = 1:64
    ch = cc(i);
    subplot(8,8,i)
    if blackoutNMS(ch,dir,mo,pos,tr) == 1
        ctl = 'fail';
%         plot((NMS(ch,:,dir,mo,pos,tr)),'r')
        plot(spktNMS(ch,:,dir,mo,pos,tr),'r')
    elseif blackoutNMS(ch,dir,mo,pos,tr) == 0
        ctl = 'pass';
%         plot((NMS(ch,:,dir,mo,pos,tr)),'b')
        plot(spktNMS(ch,:,dir,mo,pos,tr),'b')
    end
    xlim([0 8000])
%     ylim([-0.05 0.05])
    title(sprintf('%d%s',ch,ctl))
end