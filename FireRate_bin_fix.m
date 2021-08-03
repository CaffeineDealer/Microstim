function [frbl frst prCorrect fr_bl] = FireRate_bin_fix(x,blackout,Fs)

if length(size(x)) ~= 3
    disp('Check your input dimensions!')
    return
end

nch = size(x,1);
tdur = size(x,2)/2; 
nrep = size(x,3);
dur = tdur/Fs;
%% divid spike trains into X bins and create a matrix with following dimensions: nch*ndir*nmot*npos*nrep*nbin
nbin = 16;
bin = size(x,2)/nbin;
durbin = bin/Fs;

SpikeBin = zeros(nch,nrep,nbin);
for ch = 1:nch
    for t = 1:nrep
        c = 1;
        for b = 1:nbin
            SpikeBin(ch,t,b) = sum(squeeze(x(ch,c:b*bin,t)));
            c = c + bin;
        end
    end
end
SpikeBin = SpikeBin ./ durbin; % Firing Rate per each bin
%% Compute the mean of each bin across different Channel excluding the one that failed cleaning filter
for t = 1:nrep
    a = squeeze(SpikeBin(1:32,t,:));
    b = squeeze(SpikeBin(33:64,t,:));
    frmeanvP1 = mean(a(blackout(1:32)==0,:));
    frmeanvP2 = mean(b(blackout(33:64)==0,:));
    SpikeBin(1:32,t,:) = squeeze(SpikeBin(1:32,t,:)) - frmeanvP1;
    SpikeBin(33:64,t,:) = squeeze(SpikeBin(33:64,t,:)) - frmeanvP2;
end
% This would generate bunch of negative firing rates replace them with zero
SpikeBin(SpikeBin(:) < 0) = 0;
SpikeBin = permute(SpikeBin,[1 3 2]);
%%
xbl = SpikeBin(:,1:(nbin/2),:);
xst = SpikeBin(:,(nbin/2)+1:end,:);

fr_bl = zeros(nch,1);
frst = zeros(nch,1);

for ch = 1:nch
    wbl = squeeze(xbl(ch,:,:));
    wst = squeeze(xst(ch,:,:));
    bo = squeeze(blackout(ch,:));
    wbl = wbl(:,bo==0);
    wst = wst(:,bo==0);
    prCorrect(ch) = 100 - (sum(bo) * 25);
    if ~isempty(wbl) || ~isempty(wst)
        fr_bl(ch) = mean(wbl(:))/dur;
        frst(ch) = mean(wst(:))/dur;
    end
    wbl = [];
    wst = [];
    bo = [];
    bl = fr_bl(ch,squeeze(prCorrect(ch)~=0));
    frbl(ch,1) = sum(bl) / size(bl,2);
    bl = [];
end
