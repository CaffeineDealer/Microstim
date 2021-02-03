function [NMS MS] = loadspk(spka,spkb,dir,mot,pos,segment)

if size(spka,1) == size(spkb,1)
    Lspk = size(spka,1);
    a = spka(1:Lspk,1:dir,1:mot,1:pos,:); % recorded file a
    b = spkb(1:Lspk,1:dir,1:mot,1:pos,:); % recorded file b
else
    abdiff = abs(size(spka,1) - size(spkb,1))+1;
    Lspk = min(size(spka,1),size(spkb,1));
    switch segment
        case 'stim'
            a = spka(1:Lspk,1:dir,1:mot,1:pos,:); % recorded file a
            b = spkb(1:Lspk,1:dir,1:mot,1:pos,:); % recorded file b
        case 'bl'
            if size(spka,1) > size(spkb,1)
                a = spka(abdiff:end,1:dir,1:mot,1:pos,:); 
                b = spkb(:,1:dir,1:mot,1:pos,:); 
            elseif size(spka,1) < size(spkb,1)
                a = spka(:,1:dir,1:mot,1:pos,:); 
                b = spkb(abdiff:end,1:dir,1:mot,1:pos,:); 
            end
    end
end

NMS = zeros(Lspk,dir,mot,pos,size(a,5)); % No Micro-Stim trials
MS = zeros(Lspk,dir,mot,pos,size(a,5)); % Micro-Stim trials

NMS(:,[1,3,5,7],:,:,:) = b(:,[1,3,5,7],:,:,:);
NMS(:,[2,4,6,8],:,:,:) = a(:,[2,4,6,8],:,:,:);

MS(:,[1,3,5,7],:,:,:) = a(:,[1,3,5,7],:,:,:);
MS(:,[2,4,6,8],:,:,:) = b(:,[2,4,6,8],:,:,:);
