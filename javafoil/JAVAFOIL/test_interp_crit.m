clc; clearvars -except F*

load('data\measured_data\awa_100\cT_1D.mat')
addpath lib;

X            = cell(2, 1);
[X{1}, X{2}] = ndgrid(data.AWA, data.sheeting_angle);
V            = data.cT;
AWA          = deg2rad(99);
sheet_angle  = deg2rad(-81);
x            = [AWA, sheet_angle];
f            = @(x) x*-10;

cT = interp_criterion(X, V, x, 'linear', f);