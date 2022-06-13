%% 1D
clear

% Load struct
load('data\measured_data\awa_100\cT_1D.mat');
% load('data\measured_data\awa_pm_45\cT_1D.mat');

% Select operating point
AWA           = deg2rad(100);
sheet_angle_0 = deg2rad(-85);

[~, i] = min(abs(AWA - data.AWA));
[~, j] = min(abs(sheet_angle_0 - data.sheeting_angle));

HESS = gradient(gradient(data.cT));
fprintf("Hessian inverse estimate = %f\n", 1/HESS(i,j));