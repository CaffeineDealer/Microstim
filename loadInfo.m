function discardCell = loadInfo(pathInfo,type,option)

cd(pathInfo)
S = load(sprintf('S%s%s.mat',type,option));
N = load(sprintf('N%s%s.mat',type,option));
% switch option
%     case ''
%         names = {'ytu310a','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a'};
%     case 'MS'
%         names = {'ytu310a','ytu316a','ytu321a','ytu323a','ytu329a','ytu333a','ytu335a'};
% end
disCell = struct;
for i = 1:size(N.info,2)
    disCell(i).log(:,1) = S.info(i).ch; % Ch ID
    disCell(i).log(:,2) = S.info(i).discard; % MS
    disCell(i).log(:,3) = N.info(i).discard; % NoMS
end

discardCell = disCell;