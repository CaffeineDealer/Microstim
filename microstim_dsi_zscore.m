%% User defined values
clear all; clc
close all
path_firing = 'E:\MT_MST\SuperTuneFiringMatrix\ms_rmvd\';
path_spktrain = 'E:\MT_MST\SuperTuneSpkTrains\ms_rmvd\';
condition = 'S';
mtype = 2;
option = '';
switch option
    case ''
%         names = {'ytu310a','ytu312b','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a'};
        names = {'ytu310a','ytu312b','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a','ytu337a'};
    case 'MS'
%         names = {'ytu310a','ytu312b','ytu316a','ytu323a','ytu333a','ytu335a','ytu337a'};
        names = {'ytu312b','ytu316a','ytu323a','ytu329a','ytu335a','ytu337a'};
end
switch mtype
    case 1
        type = 'trans';
    case 2 
        type = 'spiral';
end
Fs = 10000;
info = struct;
%%
for j = 1:length(names)
    params = load(['E:\MT_MST\Plexon\RFiles\',names{1,j},'_TrialStructure.mat']);
    switch option 
        case ''
%             load(['E:\MT_MST\Microstim\Cell_Lib\MS_removed\Candidate Cells\','CLib_',names{1,j}(1:end-1),'.mat'])
            load(['E:\MT_MST\Microstim\Cell_Lib\MST_MU_MT_SU\MT\','CLib_',names{1,j}(1:end-1),'.mat'])
        case 'MS'
%             load(['E:\MT_MST\Microstim\Cell_Lib\MS_removed\MS Cells\','CLib_',names{1,j}(1:end-1),'.mat']) % MUA-cleaned
            load(['E:\MT_MST\Microstim\Cell_Lib\MST_MU_MT_SU\MST\','CLib_',names{1,j}(1:end-1),'.mat'])
    end
    chidx = CLib;
    clear CLib
    units = chidx;
    units(:,2) = 1;
    Dprime = [];
    bl = [];
    PD = [];
    ND = [];
    discard = [];
    FRp = []; FRn = [];
    f = [];
    b = [];
    pt = [];
    gch = [];
    df = [];
    chid = [];
    xr = [];
    for ci = 1:size(units,1)
        ch = units(ci,1);
        u = units(ci,2);
        firing = load([path_firing,names{1,j}(1,1:end-1),condition,num2str(ch),num2str(u),'firingMat.mat']);
        firing1 = squeeze(firing.firing(:,mtype,:));
        load([path_spktrain,names{1,j}(1,1:end-1),condition,num2str(ch),num2str(u),'spktrain.mat']);
        load([path_spktrain,names{1,j}(1,1:end-1),condition,num2str(ch),num2str(u),'spktrain_bl.mat']);
%         baseline = squeeze(sum(spktrainbl.spktrain_bl,1))*Fs/size(spktrainbl.spktrain_bl,1);
%         allstimfir = squeeze(sum(spktrain.spktrain,1))*Fs/size(spktrain.spktrain,1);
        firing_bl = load([path_firing,names{1,j}(1,1:end-1),condition,num2str(ch),num2str(u),'firingMat_bl.mat']);
        baseline = firing_bl.firing_bl;
        allstimfir = firing.firing;
%         [h,p] = ttest(baseline(:),allstimfir(:));
%         keepCriteria = (p <= 0.05)
        p = 0;
        keepCriteria = 1;
        bl = mean(baseline(:));
        [maxfir,I] = max(firing1(:));
        [dir,pos] = ind2sub(size(firing1),I);
        dirfir = squeeze(firing1(:,pos));
        f = [f, maxfir]; b = [b, bl];
        info(j).sp(ci).spktrain = squeeze(spktrain(:,dir,mtype,pos,:));
        info(j).sp(ci).spktrain_bl = squeeze(spktrain_bl(:,dir,mtype,pos,:));
        if keepCriteria == 1
            disp(sprintf('p = %d --> Pass!',p))
            if maxfir ~= 0
                PrefdirD = dirfir(dir);
                if dir <= params.file.taskDialogValues.numberOfDirections/2
                    NulldirD = dirfir(dir + (params.file.taskDialogValues.numberOfDirections/2));
                else
                    NulldirD = dirfir(dir - (params.file.taskDialogValues.numberOfDirections/2));
                end
                dprimeDots = (PrefdirD - NulldirD)/((PrefdirD - bl) + (NulldirD - bl));
                PD = [PD, dir];
                if dir <= params.file.taskDialogValues.numberOfDirections/2
                    ND = [ND, (dir+(params.file.taskDialogValues.numberOfDirections/2))];
                else
                    ND = [ND, (dir-(params.file.taskDialogValues.numberOfDirections/2))];
                end
                discard = [discard, nan];
                pt = [pt, pos];
                gch = [gch, ch];
                df = [df dirfir];
            elseif maxfir == 0
                warning('!')
                discard = [discard, ci];
                PD = [PD, nan];
                ND = [ND, nan];
                dprimeDots = nan;
                PrefdirD = nan;
                NulldirD = nan;
                bl = nan;
                pt = [pt, nan];
                df = [df,nan(8,1)];
            end
        elseif keepCriteria == 0
            disp(sprintf('p = %d --> Fail!',p))
            warning('!!')
            discard = [discard, ci];
            PD = [PD, nan];
            ND = [ND, nan];
            dprimeDots = nan;
            PrefdirD = nan;
            NulldirD = nan;
            bl = nan;
            pt = [pt, nan];
            df = [df,nan(8,1)];
        end
        Dprime = [Dprime,dprimeDots];
        FRp = [FRp, PrefdirD];
        FRn = [FRn, NulldirD];
        bl = [bl,bl];
        chid = [chid ch];
        xr(ci).firing = firing1;
    end
    info(j).discard = discard;
    info(j).df = df;
    info(j).xr = xr;
    info(j).Dir(:,1) = PD';
    info(j).Dir(:,2) = ND';
    info(j).DSI = Dprime;
    info(j).bl = bl;
    info(j).FR(:,1) = FRp;
    info(j).FR(:,2) = FRn;
    info(j).ctl(:,1) = f;
    info(j).ctl(:,2) = b;
    info(j).pt = pt';
    info(j).gch = gch';
    info(j).fname = names{1,j}(1,1:end-1);
    info(j).condition = condition;
    info(j).ch = chid';
    info(j).fid = size(units,1);
end
save(['E:\MT_MST\Microstim\PSTH\MUA_SUA\',sprintf('%s%s%s.mat',condition,type,option)],'info')
soundsc(rand(3000,1))
