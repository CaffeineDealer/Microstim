function [frbl frst prCorrect fr_bl] = FireRate_bin(x,blackout,Fs)

if length(size(x)) ~= 6
    disp('Check your input dimensions!')
    return
end

nch = size(x,1);
tdur = size(x,2)/2; 
ndir = size(x,3);
nmot = size(x,4);
npos = size(x,5);
nrep = size(x,6);
dur = tdur/Fs;
%% divid spike trains into X bins and create a matrix with following dimensions: nch*ndir*nmot*npos*nrep*nbin
nbin = 16;
bin = size(x,2)/nbin;
durbin = bin/Fs;

SpikeBin = zeros(nch,ndir,nmot,npos,nrep,nbin);
for ch = 1:nch
    for d = 1:ndir
        for m = 1:nmot
            for p = 1:npos
                for t = 1:nrep
                    c = 1;
                    for b = 1:nbin
                        SpikeBin(ch,d,m,p,t,b) = sum(squeeze(x(ch,c:b*bin,d,m,p,t)));
                        c = c + bin;
                    end
                end
            end
        end
    end
end
SpikeBin = SpikeBin ./ durbin; % Firing Rate per each bin
%% Compute the mean of each bin across different Channel excluding the one that failed cleaning filter
for d = 1:ndir
    for m = 1:nmot
        for p = 1:npos
            for t = 1:nrep
%                 for b = 1:nbin
                    a = squeeze(SpikeBin(1:32,d,m,p,t,:));                                    
                    b = squeeze(SpikeBin(33:64,d,m,p,t,:));
                    frmeanvP1 = mean(a(blackout(1:32,d,m,p,t)==0,:));
                    frmeanvP2 = mean(b(blackout(33:64,d,m,p,t)==0,:));
                    SpikeBin(1:32,d,m,p,t,:) = squeeze(SpikeBin(1:32,d,m,p,t,:)) - frmeanvP1;
                    SpikeBin(33:64,d,m,p,t,:) = squeeze(SpikeBin(33:64,d,m,p,t,:)) - frmeanvP2;
%                 end
            end
        end
    end
end
% This would generate bunch of negative firing rates replace them with zero
SpikeBin(SpikeBin(:) < 0) = 0;
SpikeBin = permute(SpikeBin,[1 6 2 3 4 5]);
%%
xbl = SpikeBin(:,1:(nbin/2),:,:,:,:);
xst = SpikeBin(:,(nbin/2)+1:end,:,:,:,:);

fr_bl = zeros(nch,ndir,nmot,npos);
frst = zeros(nch,ndir,nmot,npos);

for ch = 1:nch
    for d = 1:ndir
        for m = 1:nmot
            for p = 1:npos
                wbl = squeeze(xbl(ch,:,d,m,p,:));
                wst = squeeze(xst(ch,:,d,m,p,:));
                bo = squeeze(blackout(ch,d,m,p,:));
                wbl = wbl(:,bo==0);
                wst = wst(:,bo==0);
                prCorrect(ch,d,m,p) = 100 - (sum(bo) * 25);
                if ~isempty(wbl(bo==0)) && ~isempty(wst(bo==0))
                    fr_bl(ch,d,m,p) = mean(wbl(:))/dur; % I think I missed subtracting by length
                    frst(ch,d,m,p) = mean(wst(:))/dur;
                end
                wbl = [];
                wst = [];
                bo = [];
            end
        end
    end
    bl = fr_bl(ch,squeeze(prCorrect(ch,:,:,:)~=0));
    frbl(ch,1) = sum(bl) / size(bl,2);
    bl = [];
end
