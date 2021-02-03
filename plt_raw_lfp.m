function [b] = plt_raw_lfp(name)

load(['E:\MT_MST\Plexon\PDataMat\' sprintf('%sPData_despike.mat',name)])
b = PD.LFP;
flag = cellfun(@isempty,b);
b(flag==1) = [];
b = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),b,'UniformOutput',false));
t = 0:1/200:size(b,2)/200; t(end) = [];
figure
for i =1:32 
%     subplot(16,2,i)
    plot(t,b(i,:)); hold on 
end
figure
for i =1:32 
%     subplot(16,2,i) 
    plot(t,b(i+32,:)); hold on 
end