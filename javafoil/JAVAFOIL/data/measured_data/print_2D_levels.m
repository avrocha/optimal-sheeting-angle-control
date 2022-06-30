%% TACKING AWA
clear; 
% clc; close all;

load('awa_pm_45\cT_2D.mat');

cT = data.cT;
sheeting_angle_1 = data.sheeting_angle_1;
sheeting_angle_2 = data.sheeting_angle_2;

AWA = deg2rad(45);

[X, Y] = meshgrid(rad2deg(sheeting_angle_1), rad2deg(sheeting_angle_2));

[~, idx_awa] = min(abs(AWA - data.AWA));
cT        = cT(idx_awa, :, :);
cT_median = medfilt2(squeeze(cT));
diff_cT   = medfilt2(squeeze(cT)) - squeeze(cT);

cT_final                     = cT;
cT_final(abs(diff_cT) > 0.1) = cT_median(abs(diff_cT) > 0.1);

figure(1); clf(1); hold on;
surf(X', Y', squeeze(cT));
c = colorbar;
c.Label.String = 'cT';
xlabel('$\delta_1$', 'Interpreter', 'Latex')
ylabel('$\delta_2$', 'Interpreter', 'Latex')

figure(2); clf(2); hold on;
cT_temp = cT;
cT_temp(cT_temp < -2 | cT_temp > 2) = nan;
surf(X', Y', squeeze(cT_temp));
c = colorbar;
c.Label.String = 'cT';
xlabel('$\delta_1$', 'Interpreter', 'Latex')
ylabel('$\delta_2$', 'Interpreter', 'Latex')

figure(3); clf(3); hold on;
surf(X', Y', cT_median);
c = colorbar;
c.Label.String = 'cT';
xlabel('$\delta_1$', 'Interpreter', 'Latex')
ylabel('$\delta_2$', 'Interpreter', 'Latex')

figure(4); clf(4); hold on;
surf(X', Y', diff_cT);
c = colorbar;
c.Label.String = 'cT';
xlabel('$\delta_1$', 'Interpreter', 'Latex')
ylabel('$\delta_2$', 'Interpreter', 'Latex')

figure(5); clf(5); hold on;
surf(X', Y', squeeze(cT_final));
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
