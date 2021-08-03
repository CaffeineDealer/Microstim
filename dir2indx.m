function [indx] = dir2indx(dir)





d = [1:8;0:45:315]';
for i = 1:size(dir,1)
    indx(i,1) = d(d(:,2) == dir(i,1),1);
end