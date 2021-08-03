clear all; clear cache; close all; clc
cd E:\MT_MST\Plexon\PLXfiles
%% Merge a & b + create spikewaveforms
name = 'ytu337c.mat';
a = load(name);
Fs = a.params.Fsds;
allcomb = a.params.comb.allcomb;
Xa = {a.xsorted.xfs}';
spkt = [];
angles = sort(unique(allcomb(:,1)));
S = angles(logical([1 0 1 0 1 0 1 0]));
N = angles(~logical([1 0 1 0 1 0 1 0]));
j = 1; z = 1;
for i = 1:length(allcomb)
    if ~isempty(Xa{i})
        if any(allcomb(i,1) == S)
            spkt.MS(:,:,j) = Xa{i};
            j = j + 1;
        elseif any(allcomb(i,1) == N)
            spkt.NMS(:,:,z) = Xa{i};
            z = z + 1;
        end
    end
end
%% Creating a 3D matrix ch * npoints * repetition
NMS = spkt.NMS;
MS = spkt.MS;
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
[thrnms thrmaxnms] = findthreshold_fix(NMS,stdmin,stdmax);
[thrms thrmaxms] = findthreshold_fix(MS,stdmin,stdmax);
%% Generate spike trains
spktNMS = [];
spktMS = [];
nrept = size(NMS,3);
spktNMS = findspkt_fix(NMS,nrept,thrnms,thrmaxnms);%thrnms
spktMS = findspkt_fix(MS,size(MS,3),thrms,thrmaxms);
%% Firing Rate and Discarding blackouts
[frbl frst prCorrect fr_bl] = FireRate_bin_fix(spktNMS,blackoutNMS,Fs);
save(['E:\MT_MST\Microstim\MUA-Fix\',sprintf('%sN.mat',name)],'frst','frbl','spktNMS','blackoutNMS','prCorrect')
frbl = []; frst = []; prCorrect = []; fr_bl = [];
[frbl frst prCorrect fr_bl] = FireRate_bin_fix(spktMS,blackoutMS,Fs);
save(['E:\MT_MST\Microstim\MUA-Fix\',sprintf('%sS.mat',name)],'frst','frbl','spktMS','blackoutMS','prCorrect')
soundsc(rand(3000,1))
return
%% Plot Tuning Curves
dir = [0:45:315];
for i = 20:30
    tcplotMUA(dir,frst,1,3,3,frbl,i,prCorrect)
%     tcplotMUA(dir,frst,2,3,3,frbl,i,prCorrect)
end
return
%% T-Test to find good cell
clc
for ch = 1:10
    stim = squeeze(frst(ch,:,:,:));
    baseline = squeeze(fr_bl(ch,:,:,:));
    [h,p] = ttest(baseline(:),stim(:));
    if p <= 0.05
        disp('Yo! You got a Cell!')
    elseif p > 0.05
        warning('Not a Cell, Life Sucks!')
    end
end
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
