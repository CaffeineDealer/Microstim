function [NMS MS] = MS_NMS(fa,fb,dir,mot,pos)

a = fa(1:dir,1:mot,1:pos); % recorded file a
b = fb(1:dir,1:mot,1:pos); % recorded file b

NMS = zeros(dir,mot,pos); % No Micro-Stim trials
MS = zeros(dir,mot,pos); % Micro-Stim trials

NMS([1,3,5,7],:,:) = b([1,3,5,7],:,:);
NMS([2,4,6,8],:,:) = a([2,4,6,8],:,:);

MS([1,3,5,7],:,:) = a([1,3,5,7],:,:);
MS([2,4,6,8],:,:) = b([2,4,6,8],:,:);

