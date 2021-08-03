function [taskTrials,PD,ts,Fs,flag,x,goodtrials,info,new_zero,flg] = importEss(fname)

% Function importEss ... 
%%%% Inputs %%%%
% fname:    file name
%%%% Outputs %%%%
% new_zero:   offset time between MonkeyLab & Plexon

% Written by Yavar Korkian on Dec.13.2020
% Last Updated by Yavar Korkian on Dec.17.2020





load(['E:\MT_MST\Plexon\TrialStructs\' sprintf('%s_TrialStructure.mat',fname)])
cd E:\MT_MST\Plexon\PLXfiles
[~, evts3, ~] = plx_event_ts(sprintf('%s-01.plx',fname),3);
[~, evts4, ~] = plx_event_ts(sprintf('%s-01.plx',fname),4);
[~, evts6, ~] = plx_event_ts(sprintf('%s-01.plx',fname),6);
[~, evts9, ~] = plx_event_ts(sprintf('%s-01.plx',fname),9);

[~, ~, ts, ~, ~] = plx_ad_v(sprintf('%s-01.plx',fname),1);
% Modifications on May 31 for fixation only condition
% load(['E:\MT_MST\Plexon\PDataMat\' sprintf('%sPData_despike.mat',fname)]);

% Fs = PD.LFP_adfreq;
% x = PD.LFP;
Fs = 200;
x = [];
flag = [];
PD = [];
% flag = double(cellfun(@isempty,x));
% % x(flag==1) = [];
% flag = find(flag==1);
% for i = 1:size(flag,1)
%     x{flag(i)} = zeros(size(x{1},1),1);
% end
    

% x = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),x,'UniformOutput',false));

% End of modification

% Extracting trials from MonkeyLab
[goodtrials info] = getTrialML(fname,taskTrials,evts4,evts6,evts9); 

% Finding Fixation onset for the first successful trial
[new_zero flg] = fixOnset(evts3,evts4,goodtrials,ts,Fs);  


if flg == 2
    info(:,13) = info(:,13) - ts;
end
info(:,14) = info(:,13) * Fs;
% info(:,15) = ceil(info(:,14));

