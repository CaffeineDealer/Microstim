function [Vx Vy diffV] = RFcentered(S,N,badCell,type,theta,flag)





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

MST.RF = load(['E:\MT_MST\Microstim\PSTH\MUA_V1\' sprintf('mst_RFcmass_%sMS.mat',type)]);
MST.theta = load(['E:\MT_MST\Microstim\PSTH\MUA_V1\' sprintf('mst_%sMS.mat',type)]);

%%%% Center of Mass %%%%
sep = 17;
diam = 20;
c = 1;
for i = 1:max(theta(:,7))
    mtidx = find(theta(:,7)==i);
    RFx = MST.RF.Vx(MST.theta.theta(:,4)==i,1); 
    RFy = MST.RF.Vy(MST.theta.theta(:,4)==i,1); 
    for j = 1:size(RFx,1)
        cx = RFx(j,1);
        cy = RFy(j,1);
        x = [-sep+cx;cx;sep+cx;-sep+cx;cx;sep+cx;-sep+cx;cx;sep+cx];
        y = [sep+cy;sep+cy;sep+cy;cy;cy;cy;-sep+cy;-sep+cy;-sep+cy];
        for z = 1:size(mtidx,1)
                [Vx(c,1) Vy(c,1)] = cmass(fmaxS(mtidx(z),:)',x,y);
                [Vx(c,2) Vy(c,2)] = cmass(fmaxN(mtidx(z),:)',x,y);
            cc(c,1) = i; % Session
            cc(c,2) = j; % MST in session i
            cc(c,3) = z; % MT in MST session i
            c = c + 1; 
        end
    end
end
diffV(:,1) = Vx(:,1) - Vx(:,2); % x-axis
diffV(:,2) = Vy(:,1) - Vy(:,2); % y-axis

%%%% RF center of mass %%%%
figure
scatter(Vx(:,1),Vy(:,1),'FaceColor','b'); hold on
scatter(Vx(:,2),Vy(:,2),'FaceColor','r')
title(sprintf('RF^{%s}_{cm}',type))
xlabel 'Position^o'
ylabel 'Position^o'
legend({'MicroStim','No-MicroStim'},'location','NorthWest')
%%%% RF shift %%%%
figure
for i = 1:size(diffV,1)
    quiver(Vx(i,2),Vy(i,2),diffV(i,1),diffV(i,2),'LineWidth',1,'Color','k','MaxHeadSize',.6); hold on
end
title(sprintf('RF^{%s}_{cm}',type))
xlabel 'Position^o'
ylabel 'Position^o'
legend({'Shift'},'location','NorthWest')
%%%% Net arrow calculation + Plot %%%%
if flag == 1
    figure
    for i = 1:size(diffV,1)
        quiver(0,0,diffV(i,1),diffV(i,2),'LineWidth',1,'Color','b','MaxHeadSize',.6); hold on
    end
    xx = sum(diffV(:,1)); % x-axis
    yy = sum(diffV(:,2)); % y-axis
    quiver(0,0,xx,yy,'LineWidth',1.5,'Color','k','MaxHeadSize',.6); hold on
    title(sprintf('RF^{%s}_{cm}',type))
    xlabel 'Position^o'
    ylabel 'Position^o'
end