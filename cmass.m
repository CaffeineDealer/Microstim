function [Vx Vy] = cmass(r,x,y)
% Function cmass estimates the center of mass
%%%% Inputs %%%%
% r:    firing rate
% x:    x-axis coordinates / eg. x = [-17;0;17;-17;0;17;-17;0;17]; 
% y:    y-axis coordinates / eg. y = [17;17;17;0;0;0;-17;-17;-17];
%%%% Outputs %%%%
% Vx:   center of mass location on x-axis
% Vy:   center of mass location on y-axis

% Written by Yavar Korkian on Aug.03.2020

Vx = sum(r .* x) / (sum(r));
Vy = sum(r .* y) / (sum(r));


