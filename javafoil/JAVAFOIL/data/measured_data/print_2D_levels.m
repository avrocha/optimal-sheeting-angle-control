%% TACKING AWA
clear; 
% clc; close all;

load('awa_pm_45\cT_2D.mat');

cT = data.cT;
sheeting_angle_1 = data.sheeting_angle_1;
sheeting_angle_2 = data.sheeting_angle_2;

AWA = deg2rad(0);

[~, idx_awa] = min(abs(AWA - data.AWA));
cT = cT(idx_awa, :, :);
cT(cT < -2 | cT > 2) = nan;
[X, Y] = meshgrid(rad2deg(sheeting_angle_1), rad2deg(sheeting_angle_2));

figure;
surf(X', Y', squeeze(cT));
c = colorbar;
c.Label.String = 'cT';
xlabel('$\delta_1$', 'Interpreter', 'Latex')
ylabel('$\delta_2$', 'Interpreter', 'Latex')

%% AWA 100
clear; 
% clc; close all;
load('awa_100\cT_2D.mat');

cT = data.cT;
sheeting_angle_1 = data.sheeting_angle_1;
sheeting_angle_2 = data.sheeting_angle_2;

AWA = deg2rad(90);

[~, idx_awa] = min(abs(AWA - data.AWA));
cT = cT(idx_awa, :, :);
cT(cT < -2 | cT > 2) = nan;
[X, Y] = meshgrid(rad2deg(sheeting_angle_1), rad2deg(sheeting_angle_2));

figure;
surf(X', Y', squeeze(cT));
c = colorbar;
c.Label.String = 'cT';
xlabel('$\delta_1$', 'Interpreter', 'Latex')
ylabel('$\delta_2$', 'Interpreter', 'Latex')
