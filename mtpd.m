function [badCell,theta] = mtpd(N,S,svpath,svflag,option,type)





% spd = horzcat(S.df); % Micro Stim
n = 1;
for i = 1:size(N,2)
    for j = 1:size(S(i).xr,2)
        spd(:,n) = S(i).xr(j).firing(:,N(i).pt(j));
        n = n + 1;
    end
end
npd = horzcat(N.df); % No Micro Stim
mfr(:,1) = mean(spd)'; % Mean firing rate @ best position for MS
mfr(:,2) = mean(npd)'; % Mean firing rate @ best position for NMS

%%%% Form a library for cell ID and File ID + change/no change in PD %%%%
ctl = vertcat(S.ch); % Ch ID = ctl(:,1)
fid = vertcat(S.fid);
ct = 0;
for i = 1:size(fid,1)
    if i == 1 
        ctl(1:fid(i,1),2) = i; % File ID = ctl(:,2)
    else
        ctl(ct+1:ct+fid(i,1),2) = i;
    end
     ct = ct + fid(i,1);
end

%%%% Vector Average %%%%
theta = nan(size(spd,2),6);
r = nan(size(spd,2),1);
x = nan(8,size(spd,2));
for i = 1:size(spd,2)
    if ~isnan(spd(1,i)) & ~isnan(npd(1,i))
        [theta(i,1),~,r(i,1)] = polarplot(spd(:,i)); % Micro Stim
        [theta(i,2),~,r(i,1)] = polarplot(npd(:,i)); % No Micro Stim
        theta(i,3) = abs(angdiff(theta(i,1),theta(i,2))); % Degree difference
        theta(i,4) = i; % cell position
        theta(i,5) = max(spd(:,i)); % Max Firing Rate for Micro Stim
        theta(i,6) = max(npd(:,i)); % Max Firing Rate for No Micro Stim
        theta(i,7) = ctl(i,2); % cell session
        x(:,i) = spd(:,i);
        if theta(i,3) > deg2rad(20) % Define a change/no change criteria
            ctl(i,3) = 1; % Cells with change in PD
            ctl(i,6) = rad2deg(theta(i,3));
        else
            ctl(i,3) = 0; % Cells with no change in PD
            ctl(i,6) = rad2deg(theta(i,3));
        end
    else
        ctl(i,3) = nan; % Bad Cells
    end
end

badCell = find(isnan(theta(:,1)));
ctl(:,4) = vertcat(S.pt); % Micro stim spatial position
ctl(:,5) = vertcat(N.pt); % No Micro stim spatial position
theta(isnan(theta(:,1)),:) = [];
r(isnan(r(:,1))) = [];
[rv, pv] = corrcoef(theta(:,1),theta(:,2));

if svflag == 1
    save([svpath,sprintf('theta_%s_%sn.mat',type(1),option)],'theta','mfr') %ms_rmvd_150300
    save([svpath,sprintf('mt_%s%sn.mat',type,option)],'theta','spd','npd')
end

% save(['E:\MT_MST\Microstim\PSTH\MUA_SUA\',sprintf('mt_%s%sn.mat',type,option)],'theta','spd','npd')
% save(['E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\',sprintf('theta_%s_%sn.mat',type(1),option)],'theta','mfr') %ms_rmvd_150300

