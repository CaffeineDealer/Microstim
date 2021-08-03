function [xy] = extractFr(mtctl,mstctl,ctl,Fr)





for i = 1:size(ctl,1)
    idx = unique([unique(mtctl(i,:)) unique(mstctl(i,:))]);
    fr = Fr(idx,:);
    fr = mean(fr(:));
    xy(i,1) = fr;
end