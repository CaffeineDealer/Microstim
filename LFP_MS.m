clear all; clear cache; clc
set.fname = 'ytu337b';
set.analysis = 'MorletWavelet'; % Options are 'MorletWavelet', 'GrangerCausality', 'CrossCorrelation'
set.CC_type = 'AllMotion'; % CrossCorrelation over either 'AllMotion' or 'SeperateMotion'
set.win = .4; %msec
set.shift = 0; %msec2
set.mode = 1;
if set.mode == 1
    set.type = 'stim';
elseif set.mode == 2
    set.type = 'fix';
else
    error('What the fuck are you doing?')
end
set.motion_type = 3;
set.num_rept = 5;
set.ang = 8;
set.grid = 9;
set.movingwin = [0.5 0.05];
set.alpha = 0.05;
set.l = [4 8 16 25 65];
set.u = [8 12 24 55 95];
params.Fs = 200;
params.tapers = [3 5];
params.pad = 0;
params.fpass = [0 90];
params.err = [2 0.05];
% set.maxlag = 0.1 * Fs;
% y = plt_raw_lfp('ytu337a');
fname = set(1).fname;
[taskTrials,PD,ts,Fs,flag,x,goodtrials,info,new_zero,flg] = importEss(fname);
%     x = x(:,new_zero:end);  % Cutting the LFPs from desired new starting point
nch = size(x,1); % # of channels after clean up
% Grouping trials
[params.comb] = trialIdx(info);
%% MUA
% get RD files
[xout params.flag params.Fsds params.thr params.thrmax] = array2matSignal(fname);
% chop into trials
xfs = zeros(nch,(set.win * params.Fsds * 2),size(info,1));
xfix = zeros(nch,(set.win * params.Fsds),size(info,1));
xstim = xfix;
info(:,15) = ceil(info(:,13) * params.Fsds);
for j = 1:size(info,1)
    xfs(:,:,j) = xout(:,info(j,15) - (set.win * params.Fsds) + 1:info(j,15) + (set.win * params.Fsds));
    xfix(:,:,j) = xout(:,info(j,15) - (set.win * params.Fsds) + 1:info(j,15));
    xstim(:,:,j) = xout(:,info(j,15) + 1:info(j,15) + (set.win * params.Fsds));
end
t = -.4:1/params.Fsds:0.4;
t(end) = [];
% Remove mean Fixation from Stim
% mxfix = mean(mean(xfix,3),2);
% xadj = xfs - mxfix;
% Extract Signals for each condition 
xsorted = [];
for i = 1:size(params.comb.comb,2)
    a = params.comb.comb(i).idx;
    for j = 1:size(a,1)
        xsorted(i).xfs(:,:,j) = xout(:,info(a(j),15) - (set.win * params.Fsds) + 1:info(a(j),15) + (set.win * params.Fsds));
        xsorted(i).xfix(:,:,j) = xout(:,info(a(j),15) - (set.win * params.Fsds) + 1:info(a(j),15));
        xsorted(i).xstim(:,:,j) = xout(:,info(a(j),15) + 1:info(a(j),15) + (set.win * params.Fsds));
    end
end
save(['E:\MT_MST\Plexon\PLXfiles\' sprintf('%sp.mat',fname)],'params','info','xsorted')
return
%% Creat NMS & MS based on file type
params.type = fname(end);
switch params.type
    case 'a'
        if size(unique(info(:,4)),1) == 2 
            NMS = [19,36;55,72;91,108;127,144];
            MS = [1,18;37,54;73,90;109,126];
        elseif size(unique(info(:,4)),1) == 3
            NMS = [28,54;82,108;136,162;190,216];
            MS = [1,27;55,81;109,135;163,189];
        end
    case 'b'
        if size(unique(info(:,4)),1) == 2
            NMS = [1,18;37,54;73,90;109,126];
            MS = [19,36;55,72;91,108;127,144];
        elseif size(unique(info(:,4)),1) == 3
            NMS = [1,27;55,81;109,135;163,189];
            MS = [28,54;82,108;136,162;190,216];
        end
end
j = 1;
xhat = [];
for z = 1:size(NMS,1)
    for i = NMS(z,1):NMS(z,2)
        xhat(j).NMS(:,:,:) = xsorted(i).xfs(:,:,:);
        j = j + 1;
    end
end
j = 1;
for z = 1:size(MS,1)
    for i = MS(z,1):MS(z,2)
        xhat(j).MS(:,:,:) = xsorted(i).xfs(:,:,:);
        j = j + 1;
    end
end
%% Power Spectral Density using pwelch function
Pfix = [];
Pstim = [];
for i = 1:nch
    [Pstim(i,:,:),F] = pwelch(squeeze(xstim(i,:,:)),10,[],128,Fs);
    [Pfix(i,:,:),F] = pwelch(squeeze(xfix(i,:,:)),10,[],128,Fs);
    Pstim(i,:,:) = Pstim(i,:,:) .* F';
    Pfix(i,:,:) = Pfix(i,:,:) .* F';
end