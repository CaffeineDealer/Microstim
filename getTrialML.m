function [goodtrials info] = getTrialML(fname,taskTrials,evts4,evts6,evts9)

% Function getTrialML obtains experimental information from TrialStruct MonkeyLab 
%%%% Inputs %%%%
% fname:    file name
%%%% Outputs %%%%
% goodtrials
% info:   trials information

% Written by Yavar Korkian on Dec.13.2020





% pathTrialStruct =  'E:\MT_MST\Plexon\TrialStructs\'; % path
% load([pathTrialStruct sprintf('%s_TrialStructure.mat',fname)]) % load trialstructure
numtrials = length(taskTrials);

goodtrials = [];
badTrials = [];
% Find good trials & remove the first one due to delay 
goodtrials = find(arrayfun(@(taskTrials) strcmp(taskTrials.endOfTrial.eot,'correct') && ...
    ~isempty(taskTrials.endOfTrial),taskTrials))';  
% Find discarded trials
badTrials = find(arrayfun(@(taskTrials) ~strcmp(taskTrials.endOfTrial.eot,'correct') && ...
    ~isempty(taskTrials.endOfTrial),taskTrials))';

m = 1;
for j = 1:length(evts4)-1
    for k = 1:length(evts6)
        if evts6(k)>evts4(j) && evts6(k)<evts4(j+1)
            gtE6(m,1) = j;
            gtE6(m,2) = double(ismember(gtE6(m,1),goodtrials));
            m = m + 1;
            break
        end
    end
end

StimOn = [];
if size(evts9,1) - size(goodtrials,1) == 1
    goodtrials(1) = [];
    gtE6(1,:) = [];
    evts6(1) = [];
    evts6(end) = [];
    StimOn = evts6(gtE6(:,2)==1);
else
    warning('Check Trial Struct')
end

% Obtain needed information from experiment
tskTrial = {taskTrials(goodtrials).trialTaskValues};
tskTrial = [tskTrial{:}];
info = [];
info(:,1) = goodtrials;
info(:,2) = [tskTrial.dotDirection]; % Direction
info(:,3) = [tskTrial.dotDirCond];
info(:,4) = [tskTrial.dotMotion]; % 1 translation, 2 spiral, 3 shear
info(:,5) = [tskTrial.dotSpeedDeg]; % Speed
info(:,6) = [tskTrial.dotSpdCond]; % Speed
info(:,7) = [tskTrial.superTuneIndex]; % Horizontal and each row
info(:,8) = [taskTrials(goodtrials).startStimulus]; % Stim onset
info(:,9) = [taskTrials(goodtrials).stopStimulus]; % Stim offset

ftime = {taskTrials(goodtrials).fixating};
ftime = [ftime{:}];
a = [ftime.time]';
info(:,11) = floor((info(:,8) - a) * 200); % Stim onset location
ctl = a(1);
a(1) = floor((ctl - taskTrials(goodtrials(1)).fixationOn) * 200);
info(:,10) = floor((a - ctl) * 200); % Fix onset location
info(1,10) = a(1);
info(:,12) = 1:size(info,1);
info(:,13) = StimOn; % Stim onset from Plexon evets6 for good trials
