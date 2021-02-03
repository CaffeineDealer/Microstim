function spbin_avg = psth_sp(stimP,bin,ch)

stimsize = size(stimP,1);
s = 1;
e = bin - 1;
for j = 1:ceil(stimsize/bin)
    if j < ceil(stimsize/bin)
        spbin_avg(j) = sum(stimP(s:s+e,ch))/0.01;
        s = s + e + 1;
    elseif j == ceil(stimsize/bin)
        spbin_avg(j) = sum(stimP(s:end,ch))/0.01;
    end
end