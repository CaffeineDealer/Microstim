function [frbl frst prCorrect] = FireRate(x,blackout,Fs)

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

xbl = x(:,1:tdur,:,:,:,:);
xst = x(:,tdur+1:end,:,:,:,:);


trials = size(x,5);
dur = tdur/Fs;

fr_bl = zeros(nch,ndir,nmot,npos);
frst = zeros(nch,ndir,nmot,npos);

for ch = 1:nch
    for d = 1:ndir
        for m = 1:nmot
            for p = 1:npos
                wbl = sum(squeeze(xbl(ch,:,d,m,p,:)));
                wst = sum(squeeze(xst(ch,:,d,m,p,:)));
                bo = squeeze(blackout(ch,d,m,p,:));
                prCorrect(ch,d,m,p) = 100 - (sum(bo) * 25);
                ybl = sum(wbl(bo==0));
                yst = sum(wst(bo==0));
                nsucctrial = size(bo(bo==0),1);
                if ~isempty(wbl(bo==0)) && ~isempty(wst(bo==0)) 
                    fr_bl(ch,d,m,p) = ybl / (nsucctrial * dur);
                    frst(ch,d,m,p) = yst / (nsucctrial * dur);
                end
                w = [];
                b = [];
            end
        end
    end
    bl = fr_bl(ch,squeeze(prCorrect(ch,:,:,:)~=0));
    frbl(ch,1) = sum(bl) / size(bl,2);
    bl = [];
end