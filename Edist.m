function dt = Edist(x,y)

p = (x(2) - x(1))^2;
q = (y(2) - y(1))^2;

dt = sqrt(p + q);