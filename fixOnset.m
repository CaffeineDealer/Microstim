function [new_zero flg] = fixOnset(evts3,evts4,goodtrials,ts,Fs)

% Function fixOnset calculates the offeset time between MonkeyLab & Plexon start time 
%%%% Inputs %%%%
% fname:    file name
%%%% Outputs %%%%
% new_zero:   offset time between MonkeyLab & Plexon

% Written by Yavar Korkian on Dec.13.2020





% Finding Fixation onset for the first successful trial
if evts3(1) < 10
    new_zero = floor((evts4(goodtrials(1)) - ts) * Fs); % version 1
    flg = 1;
elseif evts3(1) > 10
    new_zero = floor((evts4(goodtrials(1)) - ts) * Fs); % version 2
    flg = 2;
else
    error('YoYo! Define the Plexon version, meanwhile go smoke that good kush :D ')
end