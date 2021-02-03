clear all; clear cache; close all; clc
cd E:\MT_MST\Plexon\PLXfiles
%% Merge a & b + create spikewaveforms
a = load('ytu321dp.mat');
Fs = a.params.Fsds;
allcomb = a.params.comb.allcomb;
Xa = {a.xsorted.xfs}';

spkt = [];
for i = 1:size(Xa,1)
    spkt(i).NMS = Xa{i};
end

%% Covering missing trials with zeros
for i = 1:size(spkt,2)
    if size(spkt(i).NMS,3) ~= 4
        missC = 4 - size(spkt(i).NMS,3);
        spkt(i).NMS(:,:,end + missC) = zeros;
    end
end
%% Creating a 6D matrix ch * npoints * dir * motion * position * repetition
NMS = [];
NMScell = {spkt.NMS}';
i = 1;
for j = 1:8 % #direction 
    for z = 1:size(unique(allcomb(:,2)),1) % #motion
        for w = 1:size(unique(allcomb(:,3)),1) % #position
            NMS(j,z,w,:,:,:) = double(NMScell{i,1});
            i = i + 1;
        end
    end
end
NMS = permute(NMS,[4 5 1 2 3 6]);
%% remove microstim artifact
lbthr = 0.0005;
zeroSize = 0.03 * Fs; % 30 msec 
blackoutNMS = findblackout(NMS,lbthr,zeroSize);
%% Threshold for spike sorting
stdmin = 3; % min threshold for detection
stdmax = 10; % max threshold for detection
[thrnms thrmaxnms] = findthreshold(NMS,stdmin,stdmax);
%% Generate spike trains
spktNMS = [];
nrept = 4;
spktNMS = findspkt(NMS,nrept,thrnms,thrmaxnms);%thrnms
%% Firing Rate and Discarding blackouts
[frbl frst prCorrect] = FireRate(spktNMS,blackoutNMS,Fs);
%% Plot Tuning Curves
dir = [0:45:315];
for i = 33:64
    tcplotMUA(dir,frst,1,3,3,frbl,i,prCorrect)
%     tcplotMUA(dir,frst,2,3,3,frbl,i,prCorrect)
end
return
%%
ch = 51;
mo = 2;
pos = 1;
tr = 4;
figure
for i = 1:8 
    subplot(4,2,i)
    if blackoutNMS(ch,i,mo,pos,tr) == 1
        ctl = 'fail';
        plot((NMS(ch,:,i,mo,pos,tr)),'r')
%         plot(spktNMS(ch,:,i,mo,pos,tr),'r')
    elseif blackoutNMS(ch,i,mo,pos,tr) == 0
        ctl = 'pass';
        plot((NMS(ch,:,i,mo,pos,tr)),'b')
%         plot(spktNMS(ch,:,i,mo,pos,tr),'b')
    end
    title(sprintf('%s',ctl))
end
