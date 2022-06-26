%% 1D - Hessian Estimate
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

%% 2D - Get Interpolated Hessian
clear

% Load struct
% load('data\measured_data\awa_100\cT_2D.mat');
load('data\measured_data\awa_pm_45\cT_2D.mat');

% Operating point
AWA           = deg2rad(45);
sheet_angle_0 = [deg2rad(-25), deg2rad(-25)];

% Get interpolation axes
% Grid spacing
dsa1 = abs(data.sheeting_angle_1(2) - data.sheeting_angle_1(1));
dsa2 = abs(data.sheeting_angle_2(2) - data.sheeting_angle_2(1));
% Axes
sa1 = max((sheet_angle_0(1) - dsa1), data.sheeting_angle_1(1)):dsa1:min((sheet_angle_0(1) + dsa1), data.sheeting_angle_1(end));
sa2 = max((sheet_angle_0(2) - dsa2), data.sheeting_angle_2(1)):dsa2:min((sheet_angle_0(2) + dsa2), data.sheeting_angle_2(end));

% Interpolate
cT_interp = squeeze(interpn(data.AWA, data.sheeting_angle_1, data.sheeting_angle_2, data.cT, ...
                        AWA, sa1, sa2));

% Local (numerical) hessian
[gx, gy] = gradient(cT_interp);
[gxx, gxy] = gradient(gx);
[~, gyy] = gradient(gy);

% Results
hess = [gxx(2, 2), gxy(2, 2);
        gxy(2, 2), gyy(2, 2)];

inv(hess)

%% 4D - Get Interpolated Hessian (USE SHIP STRUCTURES - [WIP])
clear
addpath lib
% Load struct
load('data\measured_data\awa_100\cT_4D.mat');
% load('data\measured_data\awa_pm_45\cT_4D.mat');

% Operating point
AWA           = deg2rad(90);
sheet_angle_0 = [deg2rad(-65.7094), deg2rad(-68.8178), deg2rad(-71.3521), deg2rad(-74.8501)]';

% Get AWA neighbors

% for each neighbor
% get interpolated axes
% interpolate cT neighborhood
% if cT neighborhood is < size(3,3) -> error
% get hessian matrix for each central element

% Get sheeting angle spacing
dsa           = data.sheeting_angle(1, 2, 1) - data.sheeting_angle(1, 1, 1);
sa_axes       = zeros(4, 3);
sa_axes(1, :) = sheet_angle_0(1)-dsa:dsa:sheet_angle_0(1)+dsa;
sa_axes(2, :) = sheet_angle_0(2)-dsa:dsa:sheet_angle_0(2)+dsa;
sa_axes(3, :) = sheet_angle_0(3)-dsa:dsa:sheet_angle_0(3)+dsa;
sa_axes(4, :) = sheet_angle_0(4)-dsa:dsa:sheet_angle_0(4)+dsa;

% Get cT neighborhood
cT_neigh = zeros(length(sa_axes(1, :)), length(sa_axes(3, :)), length(sa_axes(3, :)), length(sa_axes(4, :)));

%% [WIP]
% Criterion
J = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 2), 'cT'));    
    
for i = 1:length(sa_axes(1, :))
    for j = length(sa_axes(2, :))
        for k = length(sa_axes(3, :))
            for l = length(sa_axes(4, :))
                cT_neigh(i, j, k, l) = interp_criterion_irregular();            
            end
        end
    end
end

% [WIP]
% Local (numerical) hessian
[gx, gy] = gradient(cT_interp);
[gxx, gxy] = gradient(gx);
[~, gyy] = gradient(gy);

% Results
hess = [gxx(2, 2), gxy(2, 2);
        gxy(2, 2), gyy(2, 2)];

inv(hess);