function [N S errN errS dir] = MTalignedTC(TCN,TCS,type,site)





dir = (-180:45:135);
c = 1;
for i = 1:size(TCN,2)
    for j = 1:size(TCN(i).MST,2)
        N(:,c) = TCN(i).MST(j).MTmean;
        S(:,c) = TCS(i).MST(j).MTmean;
        c = c + 1;
    end
end

ctl = [5 0;4 6;3 7;2 8;1 0];
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
N = NN; S = SS;

mN = mean(N,2);
mS = mean(S,2);

errN = (std(N,[],2)) / (sqrt(size(N,2)));
errS = (std(S,[],2)) / (sqrt(size(S,2)));


dir = (0:45:180);
figure
switch type
    case 'trans'
        col(1) = 'k';
        col(2) = 'g';
    case 'spiral'
        col(1) = 'b';
        col(2) = 'r';
end

plot(dir,mN,col(1),'LineWidth',2)
hold on
plot(dir,mS,col(2),'LineWidth',2)
legend({'NMS','MS'})

shade(dir,mN+errN,sprintf('--%s',col(1)),dir,mN-errN,sprintf('--%s',col(1)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(1))
shade(dir,mS+errS,sprintf(':%s',col(2)),dir,mS-errS,sprintf(':%s',col(2)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(2))

% errorbar(dir,mN,errN,col(1),'LineWidth',1.5)
% hold on
% errorbar(dir,mS,errS,col(2),'LineWidth',1.5)
xlim([-10 190])
title(sprintf('Mean MT TC^{%s}',type)) % MST^{Lateral}
site = 'MST^{MS site}'; % 
xlabel(sprintf('Deviation from %s PD^o',site))
ylabel 'Firing Rate (Hz)'
