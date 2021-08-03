function [dist] = pos2dist(ref)





% Function pos2dist estimates the euclidean distance of a position from a
% reference point on the screen
%%%% Inputs %%%%
% x:    x-axis coordinates / eg. x = [-17;0;17;-17;0;17;-17;0;17]; 
% y:    y-axis coordinates / eg. y = [17;17;17;0;0;0;-17;-17;-17];
% ref:  reference coordinates / eg. ref = [17 17];
%%%% Outputs %%%%
% dist:   distance from reference point

% Written by Yavar Korkian on July.16.2021


x = [-17;0;17;-17;0;17;-17;0;17]; 
y = [17;17;17;0;0;0;-17;-17;-17];



ref = [x(ref) y(ref)];
x(:,2) = ref(1);
y(:,2) = ref(2);

temp(1,:) = [0:5];
temp(2,:) = [0,17,24,34,38,48];

for i = 1:size(x,1)
    dt(i,1) = floor(Edist(x(i,:),y(i,:)));
    dist(i,1) = temp(1,temp(2,:) == dt(i,1));
end
