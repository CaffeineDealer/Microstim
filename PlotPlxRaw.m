function [Fs,xraw] = PlotPlxRaw(pathname,fname,ch)

cd(pathname)
lb = 1;
figure
[Fs,~,~,~,ad] = plx_ad_v(fname,ch);
if ad ~= -1
    time = 0:1/Fs:(size(ad,1)-1)/Fs;
    hb = max(time);
    plot(time((Fs*lb:Fs*hb)),ad(Fs*lb:Fs*hb))
    xlim([lb hb]);
end
xlabel 'Time(sec)'
disp('done!')
xraw = ad;
