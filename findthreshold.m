function [thr thrmax] = findthreshold(x,stdmin,stdmax)

if length(size(x)) ~= 6
    disp('Check your input dimensions!')
    return
end

nch = size(x,1);
tdur = size(x,2); 
ndir = size(x,3);
nmot = size(x,4);
npos = size(x,5);
nrep = size(x,6);

x = abs(x);
lx = tdur * ndir * nmot * npos * nrep;

for ch = 1:nch
    xhat = reshape(x(ch,:,:,:,:,:),lx,1);
    NoiseStd = median(xhat)/0.6745;
    thr(ch,1) = stdmin * NoiseStd; % thr for detection
    thrmax(ch,1) = stdmax * NoiseStd; % thrmax for artifact removal
    xhat = [];
end
