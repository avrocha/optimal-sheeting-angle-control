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

HESS = gradient(gradient(data.cT(i, :), data.sheeting_angle), data.sheeting_angle);
fprintf("Hessian inverse estimate = %f\n", 1/HESS(j));

%% 2D - Get Interpolated Hessian
clear

data_source = 'tacking';

switch data_source
    case 'tacking'
        load('data\measured_data\awa_pm_45\cT_2D.mat');
        AWA           = deg2rad(45);
        sheet_angle_0 = [deg2rad(-25), deg2rad(-25)];
        
    case 'awa_100'
        load('data\measured_data\awa_100\cT_2D.mat');
        AWA           = deg2rad(100);
        sheet_angle_0 = [deg2rad(-85), deg2rad(-85)];
    otherwise
        disp('Error.')
end

% Get interpolation axes
% Grid spacing
dsa1 = abs(data.sheeting_angle_1(2) - data.sheeting_angle_1(1));
dsa2 = abs(data.sheeting_angle_2(2) - data.sheeting_angle_2(1));
dawa = abs(data.AWA(2) - data.AWA(1));
% Axes
sa1 = max((sheet_angle_0(1) - dsa1), data.sheeting_angle_1(1)):dsa1:min((sheet_angle_0(1) + dsa1), data.sheeting_angle_1(end));
sa2 = max((sheet_angle_0(2) - dsa2), data.sheeting_angle_2(1)):dsa2:min((sheet_angle_0(2) + dsa2), data.sheeting_angle_2(end));

cT = zeros(size(data.cT));
for i = 1:size(data.cT,1)
    cT(i, :, :) = medfilt2(squeeze(data.cT(i, :, :)));
end

% Interpolate
cT_interp = squeeze(interpn(data.AWA, data.sheeting_angle_1, data.sheeting_angle_2, cT, ...
                        AWA, sa1, sa2));

% Local (numerical) hessian
[gx, gy] = gradient(cT_interp, sa1, sa1);
[gxx, gxy] = gradient(gx, sa1, sa1);
[~, gyy] = gradient(gy, sa1, sa1);

% Results
hess = [gxx(2, 2), gxy(2, 2);
        gxy(2, 2), gyy(2, 2)]

inv_hess = inv(hess)

%% 4D - Get Interpolated Hessian
clear
    
data_source = 'awa_100';
    
switch data_source
    case 'tacking'
        load('data\measured_data\awa_pm_45\cT_4D.mat');
%         AWA           = deg2rad(45);
%         sheet_angle_0 = [deg2rad(-15), deg2rad(-25), deg2rad(-30), deg2rad(-35)];
            AWA           = deg2rad(75);
            sheet_angle_0 = [deg2rad(-49), deg2rad(-53), deg2rad(-57), deg2rad(-60)];
        
    case 'awa_100'
        load('data\measured_data\awa_100\cT_4D.mat');
        AWA           = deg2rad(100);
        sheet_angle_0 = [deg2rad(-77), deg2rad(-78), deg2rad(-79), deg2rad(-81)];

    otherwise
        disp('Error.')
end

% Add path for 4D median filter on data (change directory to main)
addpath lib
    
% Get AWA neighbors
[~, idx_neigh] = min(abs(data.AWA - AWA));

if data.AWA(idx_neigh) > AWA
    neighbors = [idx_neigh - 1, idx_neigh];
elseif data.AWA(idx_neigh) < AWA
    neighbors = [idx_neigh, idx_neigh + 1];
else
    neighbors = idx_neigh;
end
    
% Interpolate cT neighborhood in SA
cT_interp = zeros(length(neighbors), 3, 3, 3, 3);
sa_axes   = zeros(4, 3);
dsa       = data.sheeting_angle(1, 2, 1) - data.sheeting_angle(1, 1, 1); % Common resolution to the 4 axes
i = 1;

for idx = neighbors   
    % Run median filter on data
    cT      = squeeze(data.cT(idx, :, :, :, :));
    cT_med  = medfilt4(squeeze(data.cT(idx, :, :, :, :)));
    cT_diff = abs(cT_med - cT);
    
    cT(cT_diff > 0.05) = cT_med(cT_diff > 0.05);
    
    sa_axes(1, :) = max(sheet_angle_0(1)-dsa, data.sheeting_angle(1, 1, idx)):dsa:min(sheet_angle_0(1)+dsa, data.sheeting_angle(1, end, idx));
    sa_axes(2, :) = max(sheet_angle_0(2)-dsa, data.sheeting_angle(2, 1, idx)):dsa:min(sheet_angle_0(2)+dsa, data.sheeting_angle(2, end, idx));
    sa_axes(3, :) = max(sheet_angle_0(3)-dsa, data.sheeting_angle(3, 1, idx)):dsa:min(sheet_angle_0(3)+dsa, data.sheeting_angle(3, end, idx));
    sa_axes(4, :) = max(sheet_angle_0(4)-dsa, data.sheeting_angle(4, 1, idx)):dsa:min(sheet_angle_0(4)+dsa, data.sheeting_angle(4, end, idx));
    
    [X1, X2, X3, X4]     = ndgrid(squeeze(data.sheeting_angle(1, :, idx)), squeeze(data.sheeting_angle(2, :, idx)), ...
                            squeeze(data.sheeting_angle(3, :, idx)), squeeze(data.sheeting_angle(4, :, idx)));
    [Xq1, Xq2, Xq3, Xq4] = ndgrid(sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));

    interp_res = squeeze(interpn(X1, X2, X3, X4, squeeze(data.cT(idx, :, :, :, :)), ...
                            Xq1, Xq2, Xq3, Xq4, 'cubic')); 
        
    if any(size(interp_res) < 3)
        disp('Boundary problems in interpolation')
        cT_interp(i, :, :, :, :) = NaN(3, 3, 3, 3);
    else
        cT_interp(i, :, :, :, :) = interp_res;
    end

    i = i + 1;

end

% Interpolate neighborhood in AWA
cT_interp2 = zeros(3, 3, 3, 3);
if size(cT_interp, 1) > 1
    cT_1 = squeeze(cT_interp(1, :, :, :, :));
    cT_2 = squeeze(cT_interp(2, :, :, :, :));
    for i = 1:3^4
        cT_interp2(i) = interp1(data.AWA(neighbors), [cT_1(i), cT_2(i)], AWA);
    end
else
    cT_interp2 = cT_interp;
end

% Get Hessian
[g1, g2, g3, g4]     = gradient(cT_interp2, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
[g11, g12, g13, g14] = gradient(g1, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
[g21, g22, g23, g24] = gradient(g2, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
[g31, g32, g33, g34] = gradient(g3, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
[g41, g42, g43, g44] = gradient(g4, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));

hess = [g11(2), g12(2), g13(2), g14(2);
                    g21(2), g22(2), g23(2), g24(2);
                    g31(2), g32(2), g33(2), g34(2);
                    g41(2), g42(2), g43(2), g44(2)]

inv(hess)

%% 4D Hessian Study (AWA RANGE)
clear

AWAx = [deg2rad(61:1:79), deg2rad(86:1:124)];

hess          = zeros(4, 4, length(AWAx));
inv_hess      = zeros(size(hess));
sheet_angle_0 = zeros(length(AWAx), 4);
cT_optimal    = zeros(1, length(AWAx));

for ii = 1:length(AWAx)
    if AWAx(ii) <= deg2rad(80)
        data_source = 'tacking';
    else
        data_source = 'awa_100';
    end   
    
    switch data_source
        case 'tacking'
            load('data\measured_data\awa_pm_45\cT_4D.mat');
            
        case 'awa_100'
            load('data\measured_data\awa_100\cT_4D.mat');
    
        otherwise
            disp('Error.')
    end

    % Add path for 4D median filter on data (change directory to main)
    addpath lib
    
    optima_data = load('data\optima_calc\4D_optima.mat');     
    [~, idx] = min(abs(AWAx(ii) - optima_data.data.AWA));
    AWA           = optima_data.data.AWA(idx);
    sheet_angle_0(ii, :) = optima_data.data.X(idx, :); 
    
    % Get AWA neighbors
    [~, idx_neigh] = min(abs(data.AWA - AWA));
    
    if data.AWA(idx_neigh) > AWA
        neighbors = [idx_neigh - 1, idx_neigh];
    elseif data.AWA(idx_neigh) < AWA
        neighbors = [idx_neigh, idx_neigh + 1];
    else
        neighbors = idx_neigh;
    end
    
    % Interpolate cT neighborhood in SA
    cT_interp = zeros(length(neighbors), 3, 3, 3, 3);
    sa_axes   = zeros(4, 3);
    dsa       = data.sheeting_angle(1, 2, 1) - data.sheeting_angle(1, 1, 1); % Common resolution to the 4 axes
    i = 1;
    for idx = neighbors        
        % Run median filter on data
        cT      = squeeze(data.cT(idx, :, :, :, :));
        cT_med  = medfilt4(squeeze(data.cT(idx, :, :, :, :)));
        cT_diff = abs(cT_med - cT);
        
        cT(cT_diff > 0.05) = cT_med(cT_diff > 0.05);

        sa_axes(1, :) = max(sheet_angle_0(ii, 1)-dsa, data.sheeting_angle(1, 1, idx)):dsa:min(sheet_angle_0(ii, 1)+dsa, data.sheeting_angle(1, end, idx));
        sa_axes(2, :) = max(sheet_angle_0(ii, 2)-dsa, data.sheeting_angle(2, 1, idx)):dsa:min(sheet_angle_0(ii, 2)+dsa, data.sheeting_angle(2, end, idx));
        sa_axes(3, :) = max(sheet_angle_0(ii, 3)-dsa, data.sheeting_angle(3, 1, idx)):dsa:min(sheet_angle_0(ii, 3)+dsa, data.sheeting_angle(3, end, idx));
        sa_axes(4, :) = max(sheet_angle_0(ii, 4)-dsa, data.sheeting_angle(4, 1, idx)):dsa:min(sheet_angle_0(ii, 4)+dsa, data.sheeting_angle(4, end, idx));
        
        [X1, X2, X3, X4]     = ndgrid(squeeze(data.sheeting_angle(1, :, idx)), squeeze(data.sheeting_angle(2, :, idx)), ...
                                squeeze(data.sheeting_angle(3, :, idx)), squeeze(data.sheeting_angle(4, :, idx)));
        [Xq1, Xq2, Xq3, Xq4] = ndgrid(sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
    
        interp_res = squeeze(interpn(X1, X2, X3, X4, squeeze(data.cT(idx, :, :, :, :)), ...
                                Xq1, Xq2, Xq3, Xq4, 'cubic')); 
            
        if any(size(interp_res) < 3)
            disp('Boundary problems in interpolation')
            cT_interp(i, :, :, :, :) = NaN(3, 3, 3, 3);
        else
            cT_interp(i, :, :, :, :) = interp_res;
        end
    
        i = i + 1;
    
    end
    
    % Interpolate neighborhood in AWA
    cT_interp2 = zeros(3, 3, 3, 3);
    if size(cT_interp, 1) > 1
        cT_1 = squeeze(cT_interp(1, :, :, :, :));
        cT_2 = squeeze(cT_interp(2, :, :, :, :));
        for i = 1:3^4
            cT_interp2(i) = interp1(data.AWA(neighbors), [cT_1(i), cT_2(i)], AWA);
        end
    else
        cT_interp2 = cT_interp;
    end
    
    % Get Hessian
    [g1, g2, g3, g4]     = gradient(cT_interp2, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
    [g11, g12, g13, g14] = gradient(g1, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
    [g21, g22, g23, g24] = gradient(g2, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
    [g31, g32, g33, g34] = gradient(g3, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
    [g41, g42, g43, g44] = gradient(g4, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
    
    cT_optimal(ii) = cT_interp2(2, 2, 2, 2);
    hess(:, :, ii) = [g11(2), g12(2), g13(2), g14(2);
                        g21(2), g22(2), g23(2), g24(2);
                        g31(2), g32(2), g33(2), g34(2);
                        g41(2), g42(2), g43(2), g44(2)];
    
    inv_hess(:, :, ii) = inv(hess(:, :, ii));
end

% Temp plots
fig_cnt = 1;
save    = 1;

figure(fig_cnt); clf(fig_cnt);
% Cancellation factor - Not averaged
figure(fig_cnt); clf(fig_cnt); hold on;
sgtitle('4D Criterion | Hessian (Diag)')

subplot(4, 1, 1); hold on;
plot(rad2deg(AWAx), squeeze(hess(1, 1, :)))
xlabel('AWA [º]'), ylabel('$\hat{H}_{1, 1}$', 'Interpreter', 'Latex')

subplot(4, 1, 2); hold on;
plot(rad2deg(AWAx), squeeze(hess(2, 2, :)))
xlabel('AWA [º]'), ylabel('$\hat{H}_{2, 2}$', 'Interpreter', 'Latex')

subplot(4, 1, 3); hold on;
plot(rad2deg(AWAx), squeeze(hess(3, 3, :)))
xlabel('AWA [º]'), ylabel('$\hat{H}_{3, 3}$', 'Interpreter', 'Latex')

subplot(4, 1, 4); hold on;
plot(rad2deg(AWAx), squeeze(hess(4, 4, :)))
xlabel('AWA [º]'), ylabel('$\hat{H}_{4, 4}$', 'Interpreter', 'Latex')

if save == 1
    filename = fullfile('hess_diag.png');
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

figure(fig_cnt); clf(fig_cnt); hold on;
sgtitle('4D Criterion | Hessian (Cross)')

subplot(3, 2, 1); hold on;
plot(rad2deg(AWAx), squeeze(hess(2, 1, :)))
ylim([-1.5, 2])
xlabel('AWA [º]'), ylabel('$\hat{H}_{2, 1}$', 'Interpreter', 'Latex')

subplot(3, 2, 2); hold on;
plot(rad2deg(AWAx), squeeze(hess(3, 1, :)))
ylim([-1.5, 2])
xlabel('AWA [º]'), ylabel('$\hat{H}_{3, 1}$', 'Interpreter', 'Latex')

subplot(3, 2, 3); hold on;
plot(rad2deg(AWAx), squeeze(hess(3, 2, :)))
ylim([-1.5, 2])
xlabel('AWA [º]'), ylabel('$\hat{H}_{3, 2}$', 'Interpreter', 'Latex')

subplot(3, 2, 4); hold on;
plot(rad2deg(AWAx), squeeze(hess(4, 1, :)))
ylim([-1.5, 2])
xlabel('AWA [º]'), ylabel('$\hat{H}_{4, 1}$', 'Interpreter', 'Latex')

subplot(3, 2, 5); hold on;
plot(rad2deg(AWAx), squeeze(hess(4, 2, :)))
ylim([-1.5, 2])
xlabel('AWA [º]'), ylabel('$\hat{H}_{4, 2}$', 'Interpreter', 'Latex')

subplot(3, 2, 6); hold on;
plot(rad2deg(AWAx), squeeze(hess(4, 3, :)))
ylim([-1.5, 2])
xlabel('AWA [º]'), ylabel('$\hat{H}_{4, 3}$', 'Interpreter', 'Latex')

if save == 1
    filename = fullfile('hess_cross.png');
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

figure(fig_cnt); clf(fig_cnt); hold on;
subplot(2, 1, 1); hold on;
title('Optimal Sheeting Angles')
plot(rad2deg(AWAx), rad2deg(sheet_angle_0))
legend('Sail 1', 'Sail 2', 'Sail 3', 'Sail 4')
xlabel('AWA [º]'), ylabel('$\delta [º]$', 'Interpreter', 'Latex')

subplot(2, 1, 2); hold on;
title('Optimal cT')
plot(rad2deg(AWAx), cT_optimal)
xlabel('AWA [º]'), ylabel('cT')

if save == 1
    filename = fullfile('hess_SA.png');
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;
