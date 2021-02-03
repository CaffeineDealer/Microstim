function [spikeMat,numspikes,stimLength,baseLineLength] =nev2NFile_superTune(name,chnum,spnum)
bspath='C:\research\brainstorm_db\Test_acquisition_systems\data\Saul\@raw';
if exist(['C:\research\data\saul\' name 'events.mat'],'file')==2    
    events=load(['C:\research\data\saul\' name 'events.mat']);
else
    events=Read_headers(name);
end
load(['E:\MT_MST\Plexon\RFiles\',name,'_TrialStructure.mat'])

[a,numtrials]=size(taskTrials);
 labels={events.events.label};
 idx9=find(strcmp(labels,'Event 3904')==1);
 evts9= events.events(idx9).times; % end of trial
 idx8=find(strcmp(labels,'Event 3872')==1);
 evts8= events.events(idx8).times;% start collecting data
 idx7=find(strcmp(labels,'Event 3856')==1);
 evts7= events.events(idx7).times;  % stop stimulus
 idx6=find(strcmp(labels,'Event 3848')==1);
 evts6= events.events(idx6).times;  % stimulus on
 idx4=find(strcmp(labels,'Event 3842')==1);
 evts4= events.events(idx4).times;  % fix on
 idx3=find(strcmp(labels,'Event 3841')==1);
 evts3= events.events(idx3).times;  % start trial
 
 chspts=load([bspath,name,'.nev\',name,'.nev_ums2k_spikes\times_raw_elec_raw',num2str(chnum+256), '.mat']);
 units=unique(chspts.spikes.assigns);
 thisunit=(chspts.spikes.assigns==units(spnum));
 spts=chspts.spikes.spiketimes(thisunit);
 numspikes = length(spts);

badTrials = [];
m=1;
for j=1:numtrials
    if  ~isempty(taskTrials(j).endOfTrial) && (strcmp(taskTrials(j).endOfTrial.eot,'correct'))
        goodtrialsll(m)=j;  % goodtrials11 is the thing that actually counts - the fixation/stimulus parameters were all "good"
        stimLength(m)=taskTrials(j).endTrialState-taskTrials(j).startStimulus;
        baseLineLength(m)=taskTrials(j).startStimulus-taskTrials(j).fixationOn;
        m=m+1;
    else
        badTrials=[badTrials j];
    end
end

m=1;
for j=1:length(evts6)-1
    for k=1:length(evts9)
        if evts9(k)>evts6(j) && evts9(k)<evts6(j+1)
            goodtrials6(m)=j;   % now populate the "stimulus on" good trials
            goodtrials4(m)=k;
            m=m+1;
            break
        end
    end
end

m=m-1;  % the counter starts at 1

while m~=length(goodtrialsll)                             
    disp('Different numbers of trials in Monkeylab and Plexon data');
    
    if length(goodtrialsll)<length(goodtrials6)-1   %more than 1 trial missing
        diffTrial=goodtrials4(1:length(goodtrialsll))-goodtrialsll;
        diffTrial=find(diffTrial,1,'first');
        goodtrials4=[goodtrials4(1:diffTrial-1) goodtrials4(diffTrial+1:end)];
        goodtrials6=[goodtrials6(1:diffTrial-1) goodtrials6(diffTrial+1:end)];
    else
        goodtrials4=goodtrials4(1:end-1);   % the last trial is missing
        goodtrials6=goodtrials6(1:end-1);
    end
    m=m-1;
end

% create a column for each parameter in annular dots that contains the bin
% number so for direction min=45,step=45,num=8 you have an 8 element array
% with zero corresponding to 45 and 7 corresponding to 360 and the column
% in spike mat will indicate which bin each direction is
c=1;

spikeMat = zeros(numspikes,9);

% go through all the trials --> same as length(evts9) and number of correct
% trials 

disp(['Number of good trials ' num2str(length(goodtrialsll))]);

for i = 1:m
    % find all spike times which are between fix acquired and correct
    k=find(spts>evts4(goodtrialsll(i))&spts<evts3(goodtrialsll(i)+1));  
 
    if ~isstruct(taskTrials(goodtrialsll(i)).trialTaskValues)
       disp(['missing trialTaskValues at taskTrials(' num2str(goodtrialsll(i)) ')']);
       continue;               
    end                                                    
    for j = 0:length(k)                               
        if j==0
            % fake spike at -1000 insures that trial is represented; first entry will always be -1000 (in col 1)
            spikeMat(c,1) = -1000; 
        else
            % time of spike is when the spike occurred subtrated by the
            % "stimulus on" - these times represent time after stimulus
            spikeMat(c,1) = floor((spts(k(j)) - evts6(goodtrials6(i)))*1000);
        end                                                                  
             
        spikeMat(c,2) = i;
        spikeMat(c,3) = taskTrials(goodtrialsll(i)).trialTaskValues.dotDirection;         
        %spikeMat(c,4) = taskTrials(goodtrialsll(i)).trialTaskValues.dotDirCond;
        spikeMat(c,4) = taskTrials(goodtrialsll(i)).trialTaskValues.dotMotion; %1 traslatio, 2 spiral, 3 shear
        spikeMat(c,5) = taskTrials(goodtrialsll(i)).trialTaskValues.dotSpeedDeg;
        spikeMat(c,6) = taskTrials(goodtrialsll(i)).trialTaskValues.dotSpdCond;        
        spikeMat(c,7) = taskTrials(goodtrialsll(i)).trialTaskValues.superTuneIndex; %horizontal and each row
        % grid indices are arranged left to right up to down 
        %   1   2   3
        %   4   5   6
        %   7   8   9
        spikeMat(c,8) = taskTrials(goodtrialsll(i)).trialTaskValues.patchDiaDeg;
        spikeMat(c,9) = taskTrials(goodtrialsll(i)).trialTaskValues.dotCoherence;

        c = c + 1;
    end
   %  stimLength(i)=evts9(goodtrials(i))-evts6(goodtrials6(i));
end

lastZero=find(spikeMat(:,2)==0);
if ~isempty(lastZero)
spikeMat=spikeMat(1:lastZero-1,:);
end

disp(['saved NFile to ' NFile])
save(['E:\MT_MST\Plexon\RFiles\' name num2str(chnum) num2str(spnum) 'N.mat'] ,'spikeMat','stimLength','baseLineLength')


return;

end