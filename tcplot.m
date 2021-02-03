function tcplot(dir,x,mot,r,c,type,bl,ch,keepCriteria)

figure
h = squeeze(x(:,mot,:));
for i = 1:r*c
    plot(subplot(r,c,i),dir,squeeze(h(:,i)),'Color','b','LineWidth',1)
    hold on
    plot(dir,ones(size(dir,2),1)*bl,'Color','r','LineWidth',1)
    xlim([0 max(dir)])
    switch type
        case 'Normal'
            ylim([0 max(h(:))+2])
        case 'Z-Score'
            ylim([min(h(:)-0.2) max(h(:))+0.2])
    end
end

if keepCriteria == 1
    kc = 'Pass';
elseif keepCriteria == 0
    kc = 'Fail';
end

switch mot
    case 1
        suptitle(sprintf('Translation ch%d_{%s}',ch,kc))
    case 2
        suptitle(sprintf('Spirals ch%d_{%s}',ch,kc))
end

