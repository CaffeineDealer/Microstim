function [N S errN errS dir] = MTalignedTCmotcomb(TCNt,TCSt,TCNs,TCSs)




dir = (-180:45:135);
c = 1;
for i = 1:size(TCNt,2)
    for j = 1:size(TCNt(i).MST,2)
        Nt(:,c) = TCNt(i).MST(j).MTmean;
        St(:,c) = TCSt(i).MST(j).MTmean;
        
        Ns(:,c) = TCNs(i).MST(j).MTmean;
        Ss(:,c) = TCSs(i).MST(j).MTmean;
        
        c = c + 1;
    end
end

ctl = [5 0;4 6;3 7;2 8;1 0];

N = Nt;
S = St;
n = []; s = [];
for j = 1:size(ctl,1)
    n(1,:) = N(ctl(j,1),:);
    s(1,:) = S(ctl(j,1),:);
    if ~(ctl(j,2) == 0)
        n(2,:) = N(ctl(j,2),:);
        s(2,:) = S(ctl(j,2),:);
    end
    NN(j,:) = mean(n,1);
    SS(j,:) = mean(s,1);
    n = []; s = [];
end

N = []; S = [];
Nt = NN; St = SS;

N = Ns;
S = Ss;
n = []; s = [];
for j = 1:size(ctl,1)
    n(1,:) = N(ctl(j,1),:);
    s(1,:) = S(ctl(j,1),:);
    if ~(ctl(j,2) == 0)
        n(2,:) = N(ctl(j,2),:);
        s(2,:) = S(ctl(j,2),:);
    end
    NN(j,:) = mean(n,1);
    SS(j,:) = mean(s,1);
    n = []; s = [];
end

N = []; S = [];
Ns = NN; Ss = SS;

Nm(1,:) = max(Nt);
Nm(2,:) = max(Ns);
maxN = max(Nm);
%%
Nt = Nt ./ maxN;
St = St ./ maxN;

mNt = mean(Nt,2);
mSt = mean(St,2);

errNt = (std(Nt,[],2)) / (sqrt(size(Nt,2)));
errSt = (std(St,[],2)) / (sqrt(size(St,2)));

dir = (0:45:180);
figure

type = 'trans';
col(1) = 'k';
col(2) = 'g';

plot(dir,mNt,col(1),'LineWidth',2)
hold on
plot(dir,mSt,col(2),'LineWidth',2)
legend({'NMS','MS'})

shade(dir,mNt+errNt,sprintf('--%s',col(1)),dir,mNt-errNt,sprintf('--%s',col(1)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(1))
shade(dir,mSt+errSt,sprintf(':%s',col(2)),dir,mSt-errSt,sprintf(':%s',col(2)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(2))

% errorbar(dir,mN,errN,col(1),'LineWidth',1.5)
% hold on
% errorbar(dir,mS,errS,col(2),'LineWidth',1.5)
xlim([-10 190])
title(sprintf('Mean MT TC^{%s}',type))
xlabel(sprintf('Deviation from %s PD^o',site))
ylabel 'Firing Rate (Hz)'
%%
Ns = Ns ./ maxN;
Ss = Ss ./ maxN;

mNs = mean(Ns,2);
mSs = mean(Ss,2);

errNs = (std(Ns,[],2)) / (sqrt(size(Ns,2)));
errSs = (std(Ss,[],2)) / (sqrt(size(Ss,2)));
figure

type = 'spiral';
col(1) = 'b';
col(2) = 'r';


plot(dir,mNs,col(1),'LineWidth',2)
hold on
plot(dir,mSs,col(2),'LineWidth',2)
legend({'NMS','MS'})

shade(dir,mNs+errNs,sprintf('--%s',col(1)),dir,mNs-errNs,sprintf('--%s',col(1)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(1))
shade(dir,mSs+errSs,sprintf(':%s',col(2)),dir,mSs-errSs,sprintf(':%s',col(2)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(2))

% errorbar(dir,mN,errN,col(1),'LineWidth',1.5)
% hold on
% errorbar(dir,mS,errS,col(2),'LineWidth',1.5)
xlim([-10 190])
title(sprintf('Mean MT TC^{%s}',type))
xlabel(sprintf('Deviation from %s PD^o',site))
ylabel 'Firing Rate (Hz)'
