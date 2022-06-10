%---
% Real-data implementation of ESC controller.
% [For parallel implementations - REQUIRES PARALLEL LIBRARY]
% Gradient- and Newton-based extremum seeking controller for multivariate time-variant AWA
% J(\theta) = cT(sheeting_angle)
%---
% Copyright: Alexandre Vieira da Rocha

%% Process data
clearvars -except F*; 
clc; close all;  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

fig_cnt = 1;

% Select data source
data_source = 'awa_100';

switch data_source
    case 'tacking'
        dir      = 'data\measured_data\awa_pm_45\';
        filename = [dir, 'awa_data_45.txt'];
    case 'awa_100'
        dir      = 'data\measured_data\awa_100\';
        filename = [dir, 'awa_data_100.txt'];
    otherwise
        disp('Error: Select valid data source.\n')
end

fid = fopen(filename, 'r');
fgets(fid); % Skip title
out = fscanf(fid, ['%f', ',', '%f', ',', '%f', ',', '%f'], [4, inf]);
out = out';

time    = out(:, 1);
awa     = out(:, 2);
awa_hat = out(:, 3);
heading = out(:, 4);

% Small deviation from 5Hz
fs_data = 1/(time(2) - time(1));
n       = size(out, 1);

% Edit time to obtain constant sampling frequency from the sensor
fs_data = 5; % Both datasets are approximately 5 Hz
time    = (0:1/fs_data:(1/fs_data)*(n-1))';

% Frequency analysis from data
fx        = (0:n-1)*(fs_data/n);
y         = fft(awa);
y_hat     = fft(awa_hat);
y_abs     = abs(y).^2/n;
y_hat_abs = abs(y_hat).^2/n;

figure(fig_cnt); clf(fig_cnt);
subplot(2, 1, 1); hold on;
plot(fx, y_abs);
title('Frequency spectrum of AWA raw data')
xlabel('f [Hz]')
ylabel('|awa(jw)|')

subplot(2, 1, 2); hold on;
plot(fx, y_hat_abs);
title('Frequency spectrum of AWA filtered data')
xlabel('f [Hz]')
ylabel('|awa(jw)|')
fig_cnt = fig_cnt + 1;

%% Init controller
addpath JavaFoil;  addpath Foils; addpath lib;
global ship;

% Init configs
ship   = Ship(200);
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord

R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile, x, y, dx, chord
% R2     = Rig(75,0); % pivot x,y,  
% R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile, x, y, dx, chord

ship.addRig(R1);
% ship.addRig(R2);

ship.yaw = deg2rad(0);
scale    = calc_scale();

% Method
% 'GB' for Gradient based | 'NB' for Newton based
ES_method = 'NB';
% 1 to save figures and diary or 0 to plot figures and print diary
save = 0;

% Uncomment for parallel implementation
% p = gcp('nocreate');
% if isempty(p)
%     p = parpool('local', 4);
% end

% Simulation
fs = 10; 
dt = 1/fs;
T  = ceil(time(end));
N  = length(0:dt:T);

%
% Upsample AWA with consecutive equal samples - for loop to avoid non-int
% upsampling factors
AWA = zeros(1, N);
t_sim = (0:dt:T)';
j = 0;
for i = 1:N
    if t_sim(i) >= time(j+1); j = j + 1; end
    AWA(i) = awa(j);
end

% Initial sheeting angle
switch data_source
    case 'tacking'
        sheet_angle_0 = deg2rad(-30);
    case 'awa_100'
        sheet_angle_0 = deg2rad(-70);
end

if strcmp(ES_method, 'GB')
    fprintf("Gradient-based ESC selected\n.")

    % Uncomment params below for 1D set of parameters
    f             = 0.1; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc_hp         = 0.7*f; % HPF cutoff freq
    fc_lp         = 0.7*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = 0.0750; % gain (>0 since extremum is maximum)
 
    % Uncomment params below for 2D set of parameters
%     f             = [0.2; 0.1]; % dither freq
%     A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
%     fc_hp         = 0.09; % HPF cutoff freq 
%     fc_lp         = 0.09; % LPF cutoff freq
%     lp_bool       = false; % Use LPF
%     K             = diag([0.0750, 0.0750]); % gain (>0 since extremum is maximum)
    
    % Criterion
    J = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 1), 'cT'));
    
    % Function handle
    localShip = ship;
    gb_func   = @(J) gbesc(localShip, J, dt, N, f, A, fc_hp, fc_lp, K, sheet_angle_0, AWA, lp_bool);
    
    % Parallel execution
    F_GB = parfeval(gb_func, 3, J);

elseif strcmp(ES_method, 'NB')
    fprintf("Newton-based ESC selected\n.")

%     Uncomment params below for 1D set of parameters
    f             = 0.1; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc_hp         = 0.7*f; % HPF cutoff freq
    fc_lp         = 0.7*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = 0.0025; % gain (>0 since extremum is maximum)
    wric          = 2 * pi * 0.003; % ricatti filter parameter: 0.05f (0.01f) without (with) LPF
    ric_0         = -30;

    % Uncomment params below for 2D set of parameters
%     f             = [0.1; 0.05]; % dither freq
%     A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
%     fc_hp         = 0.7*min(f); % HPF cutoff freq
%     fc_lp         = 0.7*min(f); % LPF cutoff freq
%     lp_bool       = false; % Use LPF
%     K             = diag([0.0025, 0.0025]); % gain (>0 since extremum is maximum)
%     wric          = 2 * pi * 0.003; % ricatti filter parameter: 0.005 (0.001) without (with) LPF
%     ric_0         = diag([-30, -30]);
%     
    
    % Criterion
    J = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 2), 'cT'));

    % Function handle
    localShip = ship;
    nb_func   = @(J) nbesc(localShip, J, dt, N, f, A, fc_hp, fc_lp, K, sheet_angle_0', wric, ric_0, AWA, lp_bool);    
    
    % Parallel execution
    F_GB = parfeval(nb_func, 5, J);

else
    fprintf("Wrong method. Select either \'GB\' or \'NB\'.\n")
end

%% Plots
% Get results (Uncomment line of interest)
[sheet_angle, cT, cT_grad, cT_hessian, cT_hessian_inv] = fetchOutputs(F_GB);
[sheet_angle, cT, cT_grad, cT_hessian, cT_hessian_inv] = fetchOutputs(F_NB);

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
