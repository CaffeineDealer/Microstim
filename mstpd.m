function [badCell,theta] = mstpd(N,svpath,svflag,option,type)





npd = horzcat(N.df); % No Micro Stim

%%%% Form a library for cell ID and File ID + change/no change in PD %%%%
ctl = vertcat(N.ch); % Ch ID = ctl(:,1)
fid = vertcat(N.fid);
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
theta = nan(size(npd,2),4);
r = nan(size(npd,2),1);
x = nan(8,size(npd,2));
for i = 1:size(npd,2)
    if ~isnan(npd(1,i))
        [theta(i,1),~,r(i,1)] = polarplot(npd(:,i)); % No Micro Stim
        theta(i,2) = i; % cell position
        theta(i,3) = max(npd(:,i)); % Max Firing Rate for No Micro Stim
        theta(i,4) = ctl(i,2); % cell session
    else
        ctl(i,3) = nan; % Bad Cells
    end
end

badCell = find(isnan(theta(:,1)));
ctl(:,4) = vertcat(N.pt); % No Micro stim spatial position
theta(isnan(theta(:,1)),:) = [];
r(isnan(r(:,1))) = [];

[rv, pv] = corrcoef(theta(:,1),theta(:,2));
if svflag == 1
    save([svpath,sprintf('theta_%s_%s.mat',type(1),option)],'theta') 
    save([svpath,sprintf('mst_%s%s.mat',type,option)],'theta','npd')
end

% MUA-MST-V1V2 'E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\'
% save(['E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\',sprintf('mst_%s%s.mat',type,option)],'theta','npd') % MST1_MUA_MT_SUA
