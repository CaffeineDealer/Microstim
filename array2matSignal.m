function [xout flag Fsds thr thrmax] = array2matSignal(fname)
% [x flag Fs] = array2matSignal(fname)

RDpath = 'E:\MT_MST\Plexon\RawDataMat\';
load([RDpath sprintf('%sRawData.mat',fname)])

inpt = RD;
S = inpt.LFP;
flag = double(cellfun(@isempty,S));
% S(flag==1) = [];
flag = find(flag==1);
if ~isempty(size(S{1},1))
    for i = 1:size(flag,1)
        S{flag(i)} = zeros(size(S{1},1),1);
    end
else
    disp('Check input signal size!')
    return
end

flag = find(flag==1);
x = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),S,'UniformOutput',false));
Fs = inpt.LFP_adfreq;

% remove DC offset
x = x - mean(x,2); 
% adjust for start time
StartStimT = RD.Events.StartStimT - RD.tstart;
StopStimT = RD.Events.StopStimT - RD.tstart;
% offset = floor(RD.tstart * RD.LFP_adfreq);
% if RD.tstart * RD.LFP_adfreq > 1
%     x = x(:,offset:end);
% end
% high pass filter from 300Hz to 3KHz
[B,A] = butter(2,[300/Fs 3000/Fs]);
wbFilt = [];
for i = 1:size(x,1)
    disp(sprintf('Filtering Ch %d out of %d',i,size(x,1)))
    wbFilt(i,:) = filtfilt(B,A,x(i,:));
end
% set threshold
stdmin = 3; % min threshold for detection
stdmax = 10; % max threshold for detection
NoiseStd = median(abs(wbFilt),2)/0.6745;
thr = stdmin * NoiseStd; % thr for detection
thrmax = stdmax * NoiseStd; % thrmax for artifact removal
% downsample from 10KHz to 5KHz
% wbDws = [];
dsFactor = 1;
for i = 1:size(wbFilt,1)
    pdone = ceil(i/size(wbFilt,1) * 100);
    disp(sprintf('Processing Ch %d %d percent done',i,pdone))
    wbDws(i,:) = downsample(wbFilt(i,:),dsFactor);
end
Fsds = Fs / dsFactor;
% xout = wbDws;
xout = wbFilt;