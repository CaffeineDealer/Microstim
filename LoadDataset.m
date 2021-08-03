function [N,S] = LoadDataset(fldname,param)





type = param.type;
option = param.option;
MUA_SUActl = param.mua_sua;

switch option
    case ''
        s = load([fldname,sprintf('S%s%s.mat',type,option)]);
        S = s.info;
    case 'MS'
        S = [];
end
n = load([fldname,sprintf('N%s%s.mat',type,option)]);
N = n.info;

switch option
    case ''
        names = {'ytu310a','ytu312b','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a','ytu337a'}; % MUA_SUA
        %  names = {'ytu312b','ytu316a','ytu323a','ytu329a','ytu335a','ytu337a'};
        %  names = {'ytu310a','ytu312b','ytu316a','ytu323a','ytu333a','ytu335a'};
        %  names = {'ytu312b','ytu316a','ytu323a','ytu329a','ytu335a','ytu337a'}; % ms_rmvd_150300
        % Start MUA_SUA only
        if MUA_SUActl == 1
            N(7) = []; S(7) = [];
            N(4) = []; S(4) = [];
            N(1) = []; S(1) = [];
        end
        % end
    case 'MS'
        names = {'ytu312b','ytu316a','ytu323a','ytu329a','ytu335a','ytu337a'}; % MUA_SUA & ms_rmvd_150300
        %  names = {'ytu310a','ytu312b','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a'};        
end

% cd E:\MT_MST\Microstim\Norm\MUA_SUA\
% cd E:\MT_MST\Microstim\PSTH\ms_rmvd_150300\
% cd E:\MT_MST\Microstim\PSTH\SUA_SUA\
% cd E:\MT_MST\Microstim\PSTH\MST1_MUA_MT_SUA\
% cd E:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\ % MUA_V1


% load(sprintf('S%s%s.mat',type,option))
% S = info;
% clear info
% load(sprintf('N%s%s.mat',type,option))
% N = info;
% switch option
%     case ''
% %         names = {'ytu312b','ytu316a','ytu323a','ytu329a','ytu335a','ytu337a'};
%         names = {'ytu310a','ytu312b','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a','ytu337a'}; % MUA_SUA
% %         names = {'ytu310a','ytu312b','ytu316a','ytu323a','ytu333a','ytu335a'};
% %         names = {'ytu312b','ytu316a','ytu323a','ytu329a','ytu335a','ytu337a'}; % ms_rmvd_150300
%         % Start MUA_SUA only
%         N(7) = []; S(7) = [];
%         N(4) = []; S(4) = [];
%         N(1) = []; S(1) = [];
%         % end
%     case 'MS'
% %         names = {'ytu310a','ytu312b','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a'};
%         names = {'ytu312b','ytu316a','ytu323a','ytu329a','ytu335a','ytu337a'}; % MUA_SUA & ms_rmvd_150300
% end