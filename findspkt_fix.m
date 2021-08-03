function [spkt] = findspkt_fix(x,nrept,thr,thrmax)

% Function findspkt finds spike times from input 
%%%% Inputs %%%%
% x:    input data in form of ch*nPoints*ndir*nmot*npos*nRept
% nrept: initial number of repetions specified by user
% thr:   threshold to detect spike trains
%%%% Outputs %%%%
% spkt: output spike trains in form of ch*nPoints*nRept
% Written by Yavar Korkian on Dec.23.2020





% Check input data dimension to make sure 2nd arg is nPoints & 3rd arg nRept
if size(x,2) <= size(x,3)
    disp('Change the input dimension')
    return
end

% Computation
nch = size(x,1);
tdur = size(x,2); % trial duration timepoints
spkt = zeros(nch,tdur,nrept);

for i = 1:nch
    w = squeeze(x(i,:,:));
    for j = 1:nrept
        [pks1 locs1] = findpeaks(abs(w(:,j)),'MinPeakHeight',thr(i));
        [pks2 locs2] = findpeaks(abs(w(:,j)),'MinPeakHeight',thrmax(i));
        [loc idx] = setdiff(locs1,locs2,'Stable');
        spkt(i,loc,j) = 1;
    end
    w = [];
end

