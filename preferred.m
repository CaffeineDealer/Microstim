function [pdir,pmot,ppos,npdir] = preferred(firingst)

[maxfir,I] = max(firingst(:));
[pdir,pmot,ppos] = ind2sub(size(firingst),I);
numdir = size(firingst,1);
if pdir > numdir/2
    npdir = pdir - (numdir/2);
else
    npdir = pdir + (numdir/2);
end