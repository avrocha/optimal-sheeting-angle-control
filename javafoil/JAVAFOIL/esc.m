%---
% Gradient- and Newton-based extremum seeking controller for multivariate time-variant AWA
% J(\theta) = cT(sheeting_angle)
%---
% Copyright: Alexandre Vieira da Rocha

%% Init
clear; clc; close all;  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

addpath JavaFoil;  addpath Foils; addpath lib;
global ship;
fprintf('-------------------------------------------------------------\n');

% Init configs
ship   = Ship(200);
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord

R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile, x, y, dx, chord
R2     = Rig(75,0); % pivot x,y,  
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile, x, y, dx, chord

ship.addRig(R1);
ship.addRig(R2);

ship.yaw = deg2rad(0);
scale    = calc_scale();

% Criterion
J = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship), 'cT'));

% Method
% 'GB' for Gradient based or 'NB' for Newton based
ES_method = 'GB';
% 1 to save figures and diary or 0 to plot figures and print diary
save = 1;

%% Simulation
fs = 1; % sampling frequency (Hz)
dt = 1/fs;
T  = 50;
N  = length(0:dt:T);
    
AWA = deg2rad(45) * ones(1, N);
sheet_angle_0 = [deg2rad(-35); deg2rad(-35)]; % delta(0) [nx1]
% sheet_angle_0 = deg2rad(-35); % delta(0) [nx1]

if strcmp(ES_method, 'GB')
    fprintf("Gradient-based ESC selected\n.")

    % Uncomment params below for 1D set of parameters
%     f             = 0.1; % dither freq
%     A             = deg2rad(2); % dither amplitude
%     fc_hp         = 0.09; % HPF cutoff freq
%     fc_lp         = 0.09; % LPF cutoff freq
%     lp_bool       = false; % Use LPF
%     K             = 0.0750; % gain (>0 since extremum is maximum)
 
    % Uncomment params below for 2D set of parameters
    f             = [0.2; 0.1]; % dither freq
    A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
    fc_hp         = 0.09; % HPF cutoff freq 
    fc_lp         = 0.09; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = diag([0.0750, 0.0750]); % gain (>0 since extremum is maximum)
    
    localShip = ship;
    [sheet_angle, cT, cT_grad] = gbesc(localShip, J, dt, N, f, A, fc_hp, fc_lp, K, sheet_angle_0, AWA, lp_bool);

elseif strcmp(ES_method, 'NB')
    fprintf("Newton-based ESC selected\n.")

%     Uncomment params below for 1D set of parameters
%     f             = 0.1; % dither freq
%     A             = deg2rad(2); % dither amplitude
%     fc_hp         = 0.09; % HPF cutoff freq
%     fc_lp         = 0.09; % LPF cutoff freq
%     lp_bool       = false; % Use LPF
%     K             = 0.0025; % gain (>0 since extremum is maximum)
%     wric          = 2 * pi * 0.05 * f; % ricatti filter parameter: 0.05f (0.01f) without (with) LPF
%     ric_0         = -30;

    % Uncomment params below for 2D set of parameters
    f             = [0.2; 0.1]; % dither freq
    A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
    fc_hp         = 0.09; % HPF cutoff freq
    fc_lp         = 0.09; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = diag([0.0025, 0.0025]); % gain (>0 since extremum is maximum)
    wric          = 2 * pi * 0.05 * min(f); % ricatti filter parameter: 0.05f (0.01f) without (with) LPF
    ric_0         = diag([-30, -30]);
    
    localShip = ship;
    [sheet_angle, cT, cT_grad, cT_hessian, cT_hessian_inv] = nbesc(localShip, J, dt, N, f, A, fc_hp, fc_lp, K, sheet_angle_0, wric, ric_0, AWA, lp_bool);

else
    fprintf("Wrong method. Select either \'GB\' or \'NB\'.\n")
end

%% Plots
fig_cnt = 1;
n = length(f);

if f(2) > f(1)
    dir = strcat('plots\frequency_selection\',ES_method,'_ESC\fast_at_front\');
else 
    dir = strcat('plots\frequency_selection\',ES_method,'_ESC\slow_at_front\');
end

if save == 1
    fileID = fopen(strcat(dir,'diary.txt'),'w');
elseif save == 0
    fileID = 1;
end

% Data string
if strcmp(ES_method, 'GB')
    data_str = sprintf(['Params:\n'...
                        'AWA: %f\n'...
                        'HPF: fc = %f\n'...
                        'LPF: fc = %f\n'...
                        'Dithers: A = %s, f = %s\n'...  
                        'Integrator gain (diag): K = %s\n'], rad2deg(ship.yaw), fc_hp, fc_lp, num2str(A'), num2str(f'), num2str(diag(K)'));

elseif strcmp(ES_method, 'NB')
   data_str = sprintf(['Params:\n'...
                        'AWA: %f\n'...
                        'HPF: fc = %f\n'...
                        'LPF: fc = %f\n'...
                        'Dithers: A = %s, f = %s\n'...
                        'Ricatti Filter: wric = %f\n'...
                        'Initial Hessian inverse value: %s\n'...
                        'Integrator gain (diag): K = %s\n'], rad2deg(ship.yaw), fc_hp, fc_lp, num2str(A'), num2str(f'), wric, num2str(diag(ric_0)'), num2str(diag(K)'));
end

% Print / Save params
fprintf(fileID, "----%s----\n", datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
fprintf(fileID, data_str);
if save == 1
    fclose(fileID);
end

% References [WIP]
% - Use interpolation instead of nearest neighbours
% - Get data for AWA = [-180, 180]
% - Create 4D structure for 2D case
% sheet_angle_ref = (...)
% cT_ref          = (...)

% Plots
for i = 1:n
    figure(fig_cnt); clf(fig_cnt); hold on;
    title(strcat(ES_method, '-ESC | Sheeting Angle-', num2str(i)))
    plot(0:dt:T, rad2deg(sheet_angle(i, 1:end-1)), 'b-', 'Linewidth', 2)
    % Uncomment line below to include reference lines
    % plot(0:dt:T, sheet_angle_ref, 'r--', 'Linewidth', 1)
    xlabel('t (s)'), ylabel('$\delta_s$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'delta_', num2str(i),'.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
end

figure(fig_cnt); clf(fig_cnt); hold on;
title(strcat(ES_method, '-ESC | Thrust Coeff'))
plot(0:dt:T, cT, 'b-', 'Linewidth', 2)
% Uncomment line below to include reference lines
% plot(0:dt:T,  cT_ref, 'r--', 'Linewidth', 1)
xlabel('t (s)'), ylabel('$cT$', 'Interpreter', 'Latex')
if save == 1
    filename = fullfile(strcat(dir,'cT.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

for i = 1:n
    figure(fig_cnt); clf(fig_cnt); hold on;
    title(strcat(ES_method, '-ESC | Gradient Estimate-', num2str(i)))
    plot(0:dt:T, cT_grad(i, :), 'b-', 'Linewidth', 2)
    plot(0:dt:T, movmean(cT_grad(i, :), 10), 'r--', 'Linewidth', 1.3)
    xlabel('t (s)'), ylabel('$\hat{\zeta}$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'cT_grad_', num2str(i),'.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
end

if strcmp(ES_method, 'NB')
    figure(fig_cnt); clf(fig_cnt); hold on;
    title('NB-ESC | Hessian Estimate')
    plot(0:dt:T, reshape(cT_hessian, [n^2, N]), 'Linewidth', 1.5)
    xlabel('t (s)'), ylabel('$\hat{H}$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
    
    figure(fig_cnt); clf(fig_cnt); hold on;
    title('NB-ESC | Hessian Inverse Estimate')
    plot(0:dt:T, reshape(cT_hessian_inv(:, :, 2:end), [n^2, N]), 'Linewidth', 1.5)
    xlabel('t (s)'), ylabel('$\hat{H}^{-1}$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_inv.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
end

figure(fig_cnt); clf(fig_cnt); 
title(strcat(ES_method, '-ESC | AWA'))
plot(0:dt:T, rad2deg(AWA), 'b-', 'Linewidth', 2)
xlabel('t (s)'), ylabel('deg (º)')
if save == 1
    filename = fullfile(strcat(dir,'AWA.fig'));
    saveas(figure(fig_cnt), filename);
end

% Sail configuration
figure(fig_cnt); clf(fig_cnt);
Mesh2Java(sheet_angle(:, end), ship);
fig_cnt = fig_cnt + 1;