clear all; clc
pname =  'E:\MT_MST\Plexon\PLXfiles';
names = {'ytu308a-01.plx','ytu309a-01.plx','ytu310a-01.plx','ytu312a-01.plx','ytu316a-01.plx','ytu321a-01.plx','ytu323a-01.plx','ytu329a-01.plx','ytu333a-01.plx','ytu335a-01.plx','ytu337a-01.plx'};
%%
fname = 'ytu321a-01.plx';
ch = 56;
[Fs,xraw] = PlotPlxRaw(pname,fname,ch);