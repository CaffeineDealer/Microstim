function [Vx Vy diffV] = RF(S,N,badCell,type,theta,flag)





posInf = struct;
for i = 1:size(S,2)
    for j = 1:size(S(i).xr,2)
        posInf(i).fmaxposS(j,:) = max(S(i).xr(j).firing); % Max firing rate between 8 directions @ each 9 position
        posInf(i).fmaxposN(j,:) = max(N(i).xr(j).firing);
    end
end
fmaxS = vertcat(posInf.fmaxposS);
fmaxN = vertcat(posInf.fmaxposN);
fmaxS(badCell,:) = [];
fmaxN(badCell,:) = [];
frmax(:,1) = max(fmaxS');
frmax(:,2) = max(fmaxN');
%%%% Center of Mass %%%%
sepration = 17;
diam = 20;
x = [-17;0;17;-17;0;17;-17;0;17];
y = [17;17;17;0;0;0;-17;-17;-17];
for i = 1:size(fmaxS,1)
    [Vx(i,1) Vy(i,1)] = cmass(fmaxS(i,:)',x,y);
    [Vx(i,2) Vy(i,2)] = cmass(fmaxN(i,:)',x,y);
end
diffV(:,1) = Vx(:,1) - Vx(:,2); % x-axis
diffV(:,2) = Vy(:,1) - Vy(:,2); % y-axis
%%%% RF center of mass %%%%
figure
% cmap = [1 0 1;0 1 1;0 1 0;0 0 1;0 0 0;1 0 0];
% for i = 1:size(Vx,1)
%     col = cmap(ctl(i,2),:);
%     scatter(Vx(i,1),Vy(i,1),40,col,'Filled'); hold on
%     scatter(Vx(i,2),Vy(i,2),40,col,'Filled')
% end
scatter(Vx(:,1),Vy(:,1),'FaceColor','b'); hold on
scatter(Vx(:,2),Vy(:,2),'FaceColor','r')
line([0 0],[min(y) max(y)],'LineWidth',0.5,'Color','k')
line([min(x) max(Vx(:))],[0 0],'LineWidth',0.5,'Color','k')
xlim([min(x) max(Vx(:))])
line([min(x) max(x)],[0 0],'LineWidth',0.5,'Color','k')
xlim([min(x) max(x)])
ylim([min(y) max(y)])
title(sprintf('RF^{%s}_{cm}',type))
xlabel 'Position^o'
ylabel 'Position^o'
% legend({'MicroStim','No-MicroStim'},'location','NorthWest')
%%%% RF shift %%%%
figure
for i = 1:size(diffV,1)
    quiver(Vx(i,2),Vy(i,2),diffV(i,1),diffV(i,2),'LineWidth',1,'Color','k','MaxHeadSize',.6); hold on
end
line([0 0],[min(y) max(y)],'LineWidth',0.1,'Color','g')
% line([min(x) max(Vx(:))],[0 0],'LineWidth',0.1,'Color','g')
% xlim([min(x) max(Vx(:))])
line([min(x) max(x)],[0 0],'LineWidth',0.1,'Color','g')
xlim([min(x) max(x)])
ylim([min(y) max(y)])
title(sprintf('RF^{%s}_{cm}',type))
xlabel 'Position^o'
ylabel 'Position^o'
legend({'Shift'},'location','NorthWest')
%%%% Net arrow calculation + Plot %%%%
% load('pdDiff.mat')
if flag == 1
    figure
    for i = 1:size(diffV,1)
        [qs(i,1) cmap(i,:)] = qsign(theta(i,5),theta(i,6),'scatter'); % MT
        %     [qs(i,1) cmap(i,:)] = qsign(pdDiff(i,2),pdDiff(i,1),'scatter'); % MT vs MST Trans
        %     [qs(i,1) cmap(i,:)] = qsign(pdDiff(i,4),pdDiff(i,3),'scatter'); % MT vs MST Spiral
        quiver(0,0,diffV(i,1),diffV(i,2),'LineWidth',1,'Color',cmap(i,:),'MaxHeadSize',.6); hold on
    end
    xx = sum(diffV(:,1)); % x-axis
    yy = sum(diffV(:,2)); % y-axis
    quiver(0,0,xx,yy,'LineWidth',1.5,'Color','k','MaxHeadSize',.6); hold on
    line([0 0],[min(y) max(y)],'LineWidth',0.1,'Color','y')
    line([min(x) max(x)],[0 0],'LineWidth',0.1,'Color','y')
    % xlim([-abs(xx) abs(xx)])
    % ylim([-abs(yy) abs(yy)])
    xlim([-15 15])
    ylim([-15 15])
    % ylim([-abs(yy) abs(yy)])
    title(sprintf('RF^{%s}_{cm}',type))
    xlabel 'Position^o'
    ylabel 'Position^o'
end