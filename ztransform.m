function z = ztransform(x,dim,lb,ub)

y = x(:,lb:ub);

mu = mean(y,dim);
sigma = std(y,[],dim);


z = (x - mu) ./ sigma;