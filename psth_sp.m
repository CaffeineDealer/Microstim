function [spbin] = psth_sp(x,bin,Fs)
% psth_sp computes psth for input signal x

% Inputs:
% x: spike train of size (duration * repetition)
% bin: bin size of interest to compute PSTH
% Fs: sampling frequency in Hz

% Outputs:
% spbin: Produces output binned spike train of size 

% Written by Yavar Korkian on May.2020
% Last Modification: Feb.2021


% To make sure input is dur by rep
if size(x,1) < size(x,2)
    disp('Yo! Check your damn input dimensions')
    return
end


Inputsize = size(x,1);
nrep = size(x,2);
rate = bin / Fs;
s = 1;
e = bin - 1;

spbin = [];
for j = 1:ceil(Inputsize/bin)
    if j < ceil(Inputsize/bin)
        spbin(j,:) = sum(x(s:s+e,:));
        s = s + e + 1;
    elseif j == ceil(Inputsize/bin)
        spbin(j,:) = sum(x(s:end,:));
    end
end
% Average over trials
spbin = sum(spbin,2) / (nrep * rate);

