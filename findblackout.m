function [blackout] = findblackout(x,lbthr,zeroSize)


x = abs(x);
x(x <= lbthr) = 0; 
ctl = zeros(1,zeroSize);

blackout = [];
for i = 1:size(x,1)
    for j = 1:size(x,3)
        for w = 1:size(x,4)
            for z = 1:size(x,5)
                for k = 1:size(x,6)
                    out = [];
                    out = strfind(x(i,:,j,w,z,k),ctl);
                    if isempty(out)
                        blackout(i,j,w,z,k) = 0;
                    elseif ~isempty(out)
                        blackout(i,j,w,z,k) = 1;
                    end
                end
            end
        end
    end
end

                    