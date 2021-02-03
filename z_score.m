function z = z_score(x,nmot)

% x: input matrix to be normalized
% z: z-scores

x = x(:,1:nmot,:);

mu = mean(x(:));
sigma = std(x(:));

z = (x - mu) / sigma;