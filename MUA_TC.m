%%
clear all; clc
close all
pathf = 'E:\MT_MST\Microstim\MST-MUA\';
name = 'ytu335N.mat';
load([pathf,name]);

dir = [0:45:315];
for i = 1:32
    tcplotMUA(dir,frst,1,3,3,frbl,i,prCorrect)
    tcplotMUA(dir,frst,2,3,3,frbl,i,prCorrect)
end
