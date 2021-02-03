clear all; clc
path = 'E:\MT_MST\SuperTuneSpkTrains\';
names = {'ytu308S','ytu309S','ytu310S','ytu316S','ytu321S','ytu329S','ytu333S','ytu335S'};
u = 1;
mtype = 1;
bin = 100;
Fs = 10000;
PDcond = 'N'; % Y or N
for j = 1:size(names,2)
    name = names{j}; 
    condition = 'N'; % Microstim: S or No-Microstim: N
    load(['E:\MT_MST\Microstim\PSTH\',name(1:end-1),'Chunum',condition,PDcond,'.mat'])
    load(['E:\MT_MST\Microstim\PSTH\',name(1:end-1),'Pos',condition,PDcond,'.mat'])
    chu = chunum(chunum(:,2) == 1);
    for range = 1:length(chu)
        [stim stim_bl sz_spk sz_spk_bl] = spkt(path,name,condition,chu(range),u,mtype,pos(range));
        for i = 1:size(stim,2)
            sk(i,:) = psth_sp(stim,bin,i);
            sk_bl(i,:) = psth_sp(stim_bl,bin,i);
        end
        spike{j}(range,:) = mean(sk);
        spike_bl{j}(range,:) = mean(sk_bl);
        clear sk sk_bl
    end
end

for i = 1:9
    if size(spike_bl{i},2) > 41
        spike_bl{i}(:,1:size(spike_bl{i},2)-41,:) = [];
    end
end
sp = vertcat(spike{:});
sp_bl = vertcat(spike_bl{:});
sp_m = zeros(size(sp,1),82);
for i = 1:size(sp,1)
    sp_m(i,1:41) = (squeeze(mean(sp_bl(i,:,:),3)));
    sp_m(i,42:end) = (squeeze(mean(sp(i,:,:),3)));
end
save(sprintf('%s%s',condition,PDcond),'sp_m')
return
%%
load('SY.mat'); SY = sp_m;
load('SN.mat'); SN = sp_m;
load('NN.mat'); NN = sp_m;
load('NY.mat'); NY = sp_m;
%% Time-points prep
Fs = 10000; bin = 100;
figure
time1 = -4051/Fs:bin/Fs:0;
time2 = 0:bin/Fs:4051/Fs;
time = [time1 time2];
STime = find(time == 0.15);
Ston = find(time == 0);
%% PSTH for all Cells combined
S = zeros(113,82); N = zeros(113,82);
S(1:68,:) = SY;
S(69:end,:) = SN;
N(1:68,:) = NY;
N(69:end,:) = NN;

figure
plot(time,mean(S),'Color','r','LineWidth',2)                                    
hold on
plot(time,mean(N),'Color','k','LineWidth',2)    
line([time(Ston) time(Ston)],[0 max(mean(S))],'LineWidth',2,'Color','g')
line([time(STime) time(STime)],[0 max(mean(S))],'LineWidth',2,'Color','b')
xlim([-0.4 0.4])
ylim([0 max(mean(S))+1])
xlabel 'Time (ms)'
ylabel 'Mean Spike Rate (Hz)'
legend({'MicroStim','No-MicroStim'})
%% PSTH per session
b = [1,2;3,14;15,18;19,29;30,39;40,50;51,57;58,63;64,68;];
a = [1,3;4,10;11,16;17,21;22,27;28,35;36,36;37,41;42,45];
for i = 1:9
    stim{i} = vertcat(SY(b(i,1):b(i,2),:),SN(a(i,1):a(i,2),:));
    nstim{i} = vertcat(NY(b(i,1):b(i,2),:),NN(a(i,1):a(i,2),:));
end

figure
for i = 1:9
    subplot(3,3,i)
    plot(time,mean(stim{i}),'Color','r','LineWidth',1)
    hold on
    plot(time,mean(nstim{i}),'Color','k','LineWidth',1)
    xlim([-0.4 0.4])
    title(sprintf('%d Cells',size(stim{i},1)))
end
legend({'Microstim','No-Microstim'})
xlabel 'Time (ms)'
ylabel 'Mean Spike Rate (Hz)'




