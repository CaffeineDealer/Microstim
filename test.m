name = 'ytu309b';
for i = 1:64
    if exist(sprintf('%s%d1firingMat.mat',name,i)) > 0
        stimfiring = load([name,num2str(i),'1','firingMat']);
        firingst = stimfiring.firing;
        ss(i,1) = sum(isnan(firingst(:)));
    elseif exist(sprintf('%s%d1firingMat.mat',name,i)) == 0
        ss(i,1) = NaN;
    end
end
