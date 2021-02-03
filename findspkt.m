function [spkt] = findspkt(x,nrept,thr,thrmax)

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
% Adjust nRept in case of a missing trial
if size(x,6) ~= nrept
    disp('Different # of trials is detected!')
    nrept = size(x,6);
end
% Computation
nch = size(x,1);
tdur = size(x,2); % trial duration timepoints
ndir = size(x,3);
nmot = size(x,4);
npos = size(x,5);
spkt = zeros(nch,tdur,ndir,nmot,npos,nrept);

for i = 1:nch
    for d = 1:ndir
        for m = 1:nmot
            for p = 1:npos
                w = squeeze(x(i,:,d,m,p,:));
                for j = 1:nrept
                    [pks1 locs1] = findpeaks(abs(w(:,j)),'MinPeakHeight',thr(i));
                    [pks2 locs2] = findpeaks(abs(w(:,j)),'MinPeakHeight',thrmax(i));
                    [loc idx] = setdiff(locs1,locs2,'Stable');
                    spkt(i,loc,d,m,p,j) = 1;
%                     spkt(i,find(abs(w(:,j)) >= thr(i) & abs(w(:,j)) < thrmax(i)),d,m,p,j) = 1;
                end
                w = [];
            end
        end
    end
end
    
