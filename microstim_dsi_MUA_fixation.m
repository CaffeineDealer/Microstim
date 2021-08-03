%% User defined values
clear all; clc
close all
path_firing = 'E:\MT_MST\Microstim\MUA-Fix\'; 
path_spktrain = 'E:\MT_MST\Microstim\MUA-Fix\';
condition = 'N';
option = '';
switch option
    case ''
        names = {'ytu316c','ytu329c','ytu335c','ytu337c'};
    case 'MS'
        names = {'ytu316c','ytu329c','ytu335c','ytu337c'};
end

for j = 1:length(names)
    switch option 
        case ''
            load(['E:\MT_MST\Microstim\Cell_Lib\MUA_V1vsV2\MST_V1\','CLib_',names{1,j}(1:end-1),'.mat']) 
            load([path_firing,names{1,j}(1,1:end-1),'c.mat',condition,'.mat'])
%             if condition == 'N'
%                 spktrainbl = spktNMS(:,1:4000,:);
%                 spktrain = spktNMS(:,4001:8000,:);
%             elseif condition == 'S'
%                 spktrainbl = spktMS(:,1:4000,:);
%                 spktrain = spktMS(:,4001:8000,:);
%             end
        case 'MS'
            load(['E:\MT_MST\Microstim\Cell_Lib\MUA_V1vsV2\MST_V2\','CLib_',names{1,j}(1:end-1),'.mat']) 
            load([path_firing,names{1,j}(1,1:end-1),condition,'.mat'])
%             spktrainbl = spktNMS(:,1:4000,:);
%             spktrain = spktNMS(:,4001:8000,:);
    end
    chidx = CLib;
    clear CLib
    units = chidx;
    units(:,2) = 1;
    for ci = 1:size(units,1)
        ch = units(ci,1);
        u = units(ci,2);
        x(j).frst(ci) = frst(ch);
        x(j).frbl(ci) = frbl(ch);
    end
end
save(['E:\MT_MST\Microstim\MUA-Fix\',sprintf('%s%s.mat',condition,option)],'x') % MUA-MST-V1V2
return
%%
clear all; clc
load('E:\MT_MST\Microstim\MUA-Fix\N.mat')
N = x; clear x
load('E:\MT_MST\Microstim\MUA-Fix\S.mat')
S = x; clear x
bl(:,1) = horzcat(N.frst)';
bl(:,2) = horzcat(N.frbl)';
bl(:,3) = horzcat(S.frbl)';
bl = mean(bl,2);
ms = horzcat(S.frst)';

figure
scatter(bl,ms,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor','k')
refline(1,0)
title(sprintf('Lateral MST Firing Rate NoStim (n=%d)',size(ms,1)))
xlabel 'Baseline (spk/sec)'
ylabel 'Microstim (spk/sec)'
xlim([0 25])
ylim([0 25])