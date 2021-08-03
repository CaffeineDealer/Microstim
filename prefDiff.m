function [pdDiff Fr] = prefDiff(datapath,bin)




%%%% Input Data Section %%%%
load([datapath 'theta_t_MS.mat']); trMSTt = theta; clear theta
load([datapath 'theta_s_MS.mat']); trMSTs = theta; clear theta
% load([datapath 'trMST_t_MS.mat']); trMSTt = trMST; clear theta
% load([datapath 'trMST_s_MS.mat']); trMSTs = trMST; clear theta
load([datapath 'theta_t_n.mat']); trMTt = theta; clear theta; frMTt = mfr; clear mfr
load([datapath 'theta_s_n.mat']); trMTs = theta; clear theta; frMTs = mfr; clear mfr

%%%% Computation Section %%%%
c = 1;
for i = 1:max(trMTt(:,7))
    mtidx = find(trMTt(:,7)==i);
%     mstT = trMSTt(i);
%     mstS = trMSTs(i);
    mstT = trMSTt(trMSTt(:,4)==i,1);
    mstTFr = trMSTt(trMSTt(:,4)==i,3);
    mstS = trMSTs(trMSTs(:,4)==i,1);
    mstSFr = trMSTs(trMSTs(:,4)==i,3);
    for j = 1:size(mstT,1)
        for z = 1:size(mtidx,1)
            pdDiff(c,1) = abs(angdiff(trMTt(mtidx(z),1),mstT(j))); % MST Trans NoMS vs. MT Trans MS
            pdDiff(c,2) = abs(angdiff(trMTt(mtidx(z),2),mstT(j))); % MST Trans NoMS vs. MT Trans NoMS
            pdDiff(c,3) = abs(angdiff(trMTs(mtidx(z),1),mstS(j))); % MST Spiral NoMS vs. MT Spiral MS
            pdDiff(c,4) = abs(angdiff(trMTs(mtidx(z),2),mstS(j))); % MST Spiral NoMS vs. MT Spiral NoMS
            
%             Fr(c,1) = mstTFr(j,1) - trMTt(mtidx(z),5); % MS
%             Fr(c,2) = mstTFr(j,1) - trMTt(mtidx(z),6); % NMS
%             Fr(c,3) = mstSFr(j,1) - trMTs(mtidx(z),5); % MS
%             Fr(c,4) = mstSFr(j,1) - trMTs(mtidx(z),6); % NMS

            Fr(c,1) = frMTt(mtidx(z),1); % MS
            Fr(c,2) = frMTt(mtidx(z),2); % NMS
            Fr(c,3) = frMTs(mtidx(z),1); % MS
            Fr(c,4) = frMTs(mtidx(z),2); % NMS
            
%             Fr(c,1) = trMTt(mtidx(z),5); % MS
%             Fr(c,2) = trMTt(mtidx(z),6); % NMS
%             Fr(c,3) = trMTs(mtidx(z),5); % MS
%             Fr(c,4) = trMTs(mtidx(z),6); % NMS
            Fr(c,5) = mstTFr(j,1); % MST NMS Trans
            Fr(c,6) = mstSFr(j,1); % MST NMS Spiral
            Fr(c,7) = i;
            c = c + 1;
        end
    end
end

%%%% Plot Section %%%%
[h,pt] = kstest2(pdDiff(:,1),pdDiff(:,2));
[h,ps] = kstest2(pdDiff(:,3),pdDiff(:,4));

figure
subplot(2,2,1)
F1 = cdfplot(pdDiff(:,1));
set(F1,'LineWidth',2,'Color','b')
hold on
F2 = cdfplot(pdDiff(:,2));
set(F2,'LineWidth',2,'Color','g')
legend({'MS','NMS'},'Location','East')
grid off
title 'CDF PD Diff'
xlabel ''
ylabel 'Probability'

subplot(2,2,2)
F1 = cdfplot(pdDiff(:,3));
set(F1,'LineWidth',2,'Color','k')
hold on
F2 = cdfplot(pdDiff(:,4));
set(F2,'LineWidth',2,'Color','r')
legend({'MS','NMS'},'Location','East')
grid off
title 'CDF PD Diff'
xlabel ''
ylabel ''

subplot(2,2,3)
histogram(pdDiff(:,1),bin,'FaceColor','b')
hold on
histogram(pdDiff(:,2),bin,'FaceColor','g')
title(sprintf('PD_{Diff}^{Trans} p-val = %.2f',pt))
xlabel('(MST - MT)^o_{Diff}')
ylabel('#of pairs')
legend({'MS','NMS'})
xlim([0 3.15])

subplot(2,2,4)
histogram(pdDiff(:,3),bin,'FaceColor','k')
hold on
histogram(pdDiff(:,4),bin,'FaceColor','r')
title(sprintf('PD_{Diff}^{Spiral} p-val = %.2f',ps))
legend({'MS','NMS'})
xlim([0 3.15])

figure
subplot(1,2,1)
scatter(pdDiff(:,2),pdDiff(:,1),'MarkerEdgeColor','g','MarkerFaceColor','b'); refline(1,0)
xlim([0 3.15]); ylim([0 3.15])
title('PD_{Diff}^{Trans}')
xlabel('(MST_{NoMS} - MT_{NoMS})^o')
ylabel('(MST_{NoMS} - MT_{MS})^o')

subplot(1,2,2)
scatter(pdDiff(:,4),pdDiff(:,3),'MarkerEdgeColor','r','MarkerFaceColor','k'); refline(1,0)
xlim([0 3.15]); ylim([0 3.15])
title('PD_{Diff}^{Spiral}')
xlabel('(MST_{NoMS} - MT_{NoMS})^o')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdD = zeros(size(pdDiff,1),2);
pdD(:,1) = abs(angdiff(pdDiff(:,1),pdDiff(:,2)));
pdD(:,2) = abs(angdiff(pdDiff(:,3),pdDiff(:,4)));
[h,p] = kstest2(pdD(:,1),pdD(:,2));

figure
subplot(2,1,1)
F1 = cdfplot(pdD(:,1));
set(F1,'LineWidth',2,'Color','b')
hold on
F2 = cdfplot(pdD(:,2));
set(F2,'LineWidth',2,'Color','k')
legend('Trans','Spiral','Location','East')
grid off
title 'CDF PD'
xlabel ''
ylabel 'Probability'

subplot(2,1,2)
histogram(pdD(:,1),bin,'FaceColor','b')
hold on
histogram(pdD(:,2),bin,'FaceColor','k')
xlim([0 pi])
% legend('Trans','Spiral')
title(sprintf('MT^{PD} Diff^{NoMS - MS} from MST^{PD} p-val = %.2f',pt))
xlabel('Degree (rad)')
ylabel('Paris of MT-MST')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Delta FR vs Delta PD 
frt = Fr(:,1) - Fr(:,2);
frs = Fr(:,3) - Fr(:,4);
supermin = min(min(frs),min(frt));
supermax = max(max(frs),max(frt));
offset = 0.1;

figure
subplot(1,2,1)
patch([0 pi+0.1 pi+0.1 0],[0 0 supermin-offset supermin-offset],'b','FaceAlpha',.2)
patch([0 pi+0.1 pi+0.1 0],[0 0 supermax+offset supermax+offset],'r','FaceAlpha',.2)
hold on
scatter(pdDiff(:,2),frt,'b','filled')
xlim([0 pi+0.1])
ylim([supermin-offset supermax+offset])
title 'Translation'
xlabel '{\Delta}PD^o = {\Theta}_{MST} - {\Theta}_{MST}^{NMS}'
ylabel '{\Delta}FR (spk/sec) = Fr_{MT}^{MS} - Fr_{MST}^{NMS}'

subplot(1,2,2)
patch([0 pi+0.1 pi+0.1 0],[0 0 supermin-offset supermin-offset],'b','FaceAlpha',.2)
patch([0 pi+0.1 pi+0.1 0],[0 0 supermax+offset supermax+offset],'r','FaceAlpha',.2)
hold on
scatter(pdDiff(:,4),frs,'k','filled')
xlim([0 pi+0.1])
ylim([supermin-offset supermax+offset])
title 'Spiral'
xlabel '{\Delta}PD^o = {\Theta}_{MST} - {\Theta}_{MST}^{NMS}'
