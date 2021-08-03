function [thr thrmax] = findthreshold_fix(x,stdmin,stdmax)

if length(size(x)) ~= 3
    disp('Check your input dimensions!')
    return
end

nch = size(x,1);
tdur = size(x,2); 
nrep = size(x,3); 

x = abs(x);
lx = tdur * nrep;

for ch = 1:nch
    xhat = reshape(x(ch,:,:),lx,1);
    NoiseStd = median(xhat)/0.6745;
    thr(ch,1) = stdmin * NoiseStd; % thr for detection
    thrmax(ch,1) = stdmax * NoiseStd; % thrmax for artifact removal
    xhat = [];
end