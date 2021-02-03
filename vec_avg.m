function [mean,pvalue,r] = vec_avg(V,dirs)
% rayleighz(V) returns mean angle, pvalue, and circular s.d. for vector V

% dirs = linspace(0,2*pi*(1-1/length(V)),length(V)); % matrix of : pi/4,pi/2,...,7pi/4
mean = 0.0;
pvalue = 1.0;
r = 0.0;
n = sum(V);
if n <= 0
    return;
end;
x = cos(dirs) .* V;   % cos(dirs) = cos {0,pi/4,pi/2,...,7pi/4}
y = sin(dirs) .* V;

X = sum(x)/8;        % cos mean
Y = sum(y)/8;         % sine mean

r = sqrt(X*X + Y*Y);    % circular std dev. calc.
%correction factor for large angular intervals:
d = 360/length(V);
r = r * (d*pi/360)/sin(d*pi/360);
if Y > 0 && X > 0
    mean = atan(Y/X); % mean value in radians - the preferred direction in radians.
elseif X < 0
    mean = atan(Y/X) + pi;
elseif Y < 0 && X > 0
    mean = atan(Y/X) + 2*pi;
end
% mean = atan2(Y,X);        % mean value in radians - the preferred direction in radians.
% if mean<0.0
%     mean = mean + 2*pi;
% end;
R = n * r;
pvalue = exp(sqrt(1+4*n+4*(n*n-R*R))-(1+2*n));

return;