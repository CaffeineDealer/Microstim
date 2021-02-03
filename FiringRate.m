function fr = FiringRate(x,stimst,stimend,ratio)

y = squeeze(sum(x));
y = squeeze(sum(y,4));
trials = size(x,5);
% stim = stimend - stimst;
stim = size(x,1)/10;

fr = y / (trials * stim / ratio);
% fr = squeeze(mean(fr,4);


