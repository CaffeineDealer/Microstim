function tcplotMUA(dir,x,mot,r,c,bl,ch,prCorrect)

figure
% for cc = 1:size(chh,2)
%     ch = chh(cc);
h = squeeze(x(ch,:,mot,:));
prCorrect = prCorrect ./ 25;
prc = squeeze(prCorrect(ch,:,mot,:));

for i = 1:r*c
    plot(subplot(r,c,i),dir,h(:,i),'Color','b','LineWidth',1)
%         plot(subplot(r,c,i),dir,h(:,i))
    hold on
    plot(dir,ones(size(dir,2),1)*bl(ch),'Color','r','LineWidth',1)
    xlim([0 max(dir)])
    ylim([0 max(h(:))+2])
    for j = 1:size(h,1)
        text(dir(j),floor(max(h(:)))-.1,sprintf('%d',prc(j,i)))
    end
end
% end
switch mot
    case 1
        suptitle(sprintf('Translation ch%d',ch))
    case 2
        suptitle(sprintf('Spirals ch%d',ch))
end