i = 6;
j = 2;
z = 5;
w = 6;
k = 5; 
y = 13;
%
ww = 3;
%
dirs = deg2rad(0:45:315);

figure
subplot(ww,1,1)
plot(rad2deg(dirs),trs(i).firingN(:,j),'Color','m','LineWidth',1,'LineStyle','--'); xlim([0 315])

subplot(ww,1,2)
plot(rad2deg(dirs),MT.ri_NMS_S(:,info(i).MT(z)),'Color','k','LineWidth',1,'LineStyle','--'); hold on
plot(rad2deg(dirs),MT.ri_MS_S(:,info(i).MT(z)),'Color','k','LineWidth',2); xlim([0 315])

subplot(ww,1,3)
plot(rad2deg(dirs),MT.ri_NMS_S(:,info(i).MT(w)),'Color','k','LineWidth',1,'LineStyle','--'); hold on
plot(rad2deg(dirs),MT.ri_MS_S(:,info(i).MT(w)),'Color','k','LineWidth',2); xlim([0 315])

% subplot(ww,1,4)
% plot(rad2deg(dirs),MT.ri_NMS_T(:,info(i).MT(k)),'Color','b','LineWidth',1,'LineStyle','--'); hold on
% plot(rad2deg(dirs),MT.ri_MS_T(:,info(i).MT(k)),'Color','b','LineWidth',2); xlim([0 315])

% subplot(ww,1,5)
% plot(rad2deg(dirs),MT.ri_NMS_T(:,info(i).MT(y)),'Color','b','LineWidth',1,'LineStyle','--'); hold on
% plot(rad2deg(dirs),MT.ri_MS_T(:,info(i).MT(y)),'Color','b','LineWidth',2); xlim([0 315])