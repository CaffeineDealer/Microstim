function [cor] = masterPlot(MST,MT,MTms,type)

% x


x = [];
for i = 1:size(MT,2)
    w = 1;
    for z = 1:size(MST(i).df,2)
        for j = 1:size(MT(i).df,2)
            [x(i).Fr(w,1) x(i).PD(w,1)] = max(MST(i).df(:,z)); % MST max Fr & dir
            x(i).PD(w,1) = indx2dir(x(i).PD(w,1));
            [x(i).Fr(w,2) x(i).PD(w,2)] = max(MT(i).df(:,j)); % MT max Fr & dir
            x(i).PD(w,2) = indx2dir(x(i).PD(w,2));
            w = w + 1;
            x(i).MTfr(:,j) = MT(i).df(:,j); % x(i).MT(j).firing = .firing
            x(i).MTmsfr(:,j) = MTms(i).df(:,j);
        end
    end
    x(i).MST = vertcat(MST(i).xr.firing);
end
ctl = allcomb(0:45:180,0:45:180);

for i = 1:size(x,2)
    c = 1;
    for z = 1:size(MST(i).df,2)
        for j = 1:size(MT(i).df,2)
            Frmt = x(i).MTfr(:,j); %x(i).MT(j).firing;
            Frmtms = x(i).MTmsfr(:,j); %x(i).MTms(j).firing;
            mstPD = x(i).PD(c,1);
            mtPD = x(i).PD(c,2);
            c = c + 1;
            [mstctl] = findintsec(mstPD,ctl(:,1));
            [mtctl] = findintsec(mtPD,ctl(:,2));
            mstctl(:,1) = dir2indx(mstctl(:,1));
            mstctl(:,2) = dir2indx(mstctl(:,2));
            mtctl(:,1) = dir2indx(mtctl(:,1));
            mtctl(:,2) = dir2indx(mtctl(:,2));
            x(i).dim(z,j,:) = extractFr(mtctl,mstctl,ctl,Frmt);
            x(i).dimms(z,j,:) = extractFr(mtctl,mstctl,ctl,Frmtms);
        end
    end
    x(i).xy = squeeze(mean(mean(x(i).dim,1),2));
    x(i).xyms = squeeze(mean(mean(x(i).dimms,1),2));
end

xy = mean(horzcat(x.xy),2);
xyms = mean(horzcat(x.xyms),2);

cornms = flip(reshape(xy,5,5)');
corms = flip(reshape(xyms,5,5)');

cor = corms - cornms;

figure
imagesc(cor)
% colormap jet
title(sprintf('%s',type))
xlabel 'Stim Dir from MT PD^{o}' %MST^{Lateral}
ylabel 'Stim Dir from MST PD^{o}' %^{MS Site}
xl = 0:45:180;
yl = 180:-45:0;
set(gca,'XTick',[1:1:5],'XTickLabel',xl)
set(gca,'YTick',[1:1:5],'YTickLabel',yl)
colorbar