%---
% Real-data implementation of ESC controller.
% [Interpolated criterion - Fast computations for 1D scenario]
% Gradient- and Newton-based extremum seeking controller for time-variant AWA
% J(\theta) = cT(sheeting_angle)
%---
% Copyright: Alexandre Vieira da Rocha

%% Process data
clearvars -except F*; 
clc; close all;  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

fig_cnt = 1;

% Menus
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

% Uncomment lines below to plot frequency content of AWA data
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
fprintf('-------------------------------------------------------------\n');

% Init configs
ship   = Ship(200);
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord

R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile, x, y, dx, chord

ship.addRig(R1);

ship.yaw = deg2rad(0);
scale    = calc_scale();

% Method
% 'GB' for Gradient based | 'NB' for Newton based
ES_method = 'GB';
% 1 to save figures and diary or 0 to plot figures and print diary
save = 1;

% Simulation
fs = 10; 
dt = 1/fs;
% T  = ceil(time(end));
T  = 600;
N  = length(0:dt:T);

% Upsample AWA with consecutive equal samples - for loop to avoid non-int
% upsampling factors
% AWA = zeros(1, N);
% t_sim = (0:dt:T)';
% j = 0;
% for i = 1:N
%     if t_sim(i) >= time(j+1); j = j + 1; end
%     AWA(i) = awa(j);
% end

AWA = deg2rad(100) * ones(1, N);

% Feedforward (Piecewise constant)
switch data_source
    case 'tacking'
        FF = zeros(1, N);
        [~, idx_1] = min(abs(310 - t_sim));
        [~, idx_2] = min(abs(1150 - t_sim));
        [~, idx_3] = min(abs(1560 - t_sim));
        
        FF(1:idx_1)     = deg2rad(-25);
        FF(idx_1:idx_2) = deg2rad(25);
        FF(idx_2:idx_3) = deg2rad(-25);
        FF(idx_3:end)   = deg2rad(25);

    case 'awa_100'
        FF = deg2rad(-85) * ones(1, N);
end

sheet_angle_0 = FF(1);

% Prep Interpolation
% Choose data source
switch data_source
    case 'tacking'
        load('data\measured_data\awa_pm_45\cT_1D.mat')
    case 'awa_100'
        load('data\measured_data\awa_100\cT_1D.mat')
    otherwise
        disp('Error: Select valid data source.\n')
end

X            = cell(2, 1);
[X{1}, X{2}] = ndgrid(data.AWA, data.sheeting_angle);
V            = data.cT;

% Filtering cT - {'RAW', 'EMA', 'LPF'}
cT_filter = 'RAW';
switch cT_filter
    case 'EMA' % Param = EMA 'speed'
        switch data_source
            case 'tacking'
                cT_filter_param = 0.2;   
            case 'awa_100'
                cT_filter_param = 0.1;
        end

    case 'LPF' % Param = cut-off frequency for LPF
        cT_filter_param = 0.5; 
    
    case 'RAW'
        cT_filter_param = 1;
end


if strcmp(ES_method, 'GB')
    fprintf("Gradient-based ESC selected.\n")
    
    % 1D set of parameters
    f             = 0.01; % tuning param: constant coeff
    delta         = 0.1;  % tuning param: constant coeff

    f_dither      = 20*f; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc_hp         = 3*f; % HPF cutoff freq
    fc_lp         = 5*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    
    % Data-specific tuning
    switch data_source
        case 'tacking'
            ric_0 = -0.040598;
        case 'awa_100'
            ric_0 = -0.014203;
    end
    
    K = f * delta * 10 * (-ric_0); % gain (>0 since extremum is maximum)

    % Criterion
    Jgb = @(sheeting_angle, ship) (getfield(calc_objective_mod(sheeting_angle, ship), 'cT'));
    % Interpolated criterion
    Jgb_interp = @(sheeting_angle, ship) interp_criterion(X, V, [ship.yaw, sheeting_angle'], 'linear', Jgb, ship);
    
    localShip = ship;
    [sheet_angle, cT, cT_grad, cT_hat] = gbesc(localShip, Jgb_interp, dt, N, f_dither, A, fc_hp, ...
                            fc_lp, K, sheet_angle_0, AWA, lp_bool, cT_filter, cT_filter_param, FF);   

end

if strcmp(ES_method, 'NB')
    fprintf("Newton-based ESC selected.\n")

    % 1D set of parameters
    f             = 0.01; % tuning param: constant coeff
    delta         = 0.1;  % tuning param: constant coeff

    f_dither      = 20*f; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc_hp         = 3*f; % HPF cutoff freq
    fc_lp         = 5*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = f * delta * 10; % gain (>0 since extremum is maximum)
    
    % Data-specific tuning
    switch data_source
        case 'tacking'
            ric_0 = -0.040598;
            wric  = 2 * pi * (0.01 * f * delta); % ricatti filter parameter
        case 'awa_100'
            ric_0 = -0.014203;
            wric  = 2 * pi * (10 * f * delta); % ricatti filter parameter 
    end

    % Criterion
    Jnb = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 2), 'cT'));    
    % Interpolated criterion - high resolution -> linear interpolation
    Jnb_interp = @(sheeting_angle, ship) interp_criterion(X, V, [ship.yaw, sheeting_angle'], 'linear', Jnb, ship);

    localShip = ship;
    [sheet_angle, sheet_angle_hat, cT, cT_grad, cT_hessian, cT_hessian_inv, cT_hat, hpf] = nbesc(localShip, Jnb_interp, dt, N, f_dither, A, ...
        fc_hp, fc_lp, K, sheet_angle_0, wric, ric_0, AWA, lp_bool, cT_filter, cT_filter_param, FF);    
    
end

%% Plots 
% dir = ['plots\7m_data_tacking\real\', ES_method,'_ESC\', cT_filter,'_cT\f_0_20\'];
dir = ['plots\final\AWA_100_sim\1D\', ES_method,'_ESC\'];

% Check directory
if ~exist(dir, 'dir')
    mkdir(dir)
end

n = length(f_dither);

% References
% Get reference in original grid
sheet_angle_ref = zeros(1, length(data.AWA));
cT_ref          = zeros(1, length(data.AWA ));

for i = 1:length(data.AWA)
   [cT_ref(i), idx_SA] = max(data.cT(i, :));
   sheet_angle_ref(i) = data.sheeting_angle(idx_SA);
end

% Interpolate references
sheet_angle_ref = interp1(data.AWA, sheet_angle_ref, AWA);
cT_ref          = interp1(data.AWA, cT_ref, AWA);

% Plots
for i = 1:n
    figure(fig_cnt); clf(fig_cnt); hold on;
    title(strcat(ES_method, '-ESC | Sheeting Angle-', num2str(i)))
    plot(0:dt:T, rad2deg(sheet_angle(i, 1:end-1)), 'b-', 'Linewidth', 2)
    % Uncomment line below to include reference lines
    plot(0:dt:T, rad2deg(sheet_angle_ref), 'r--', 'Linewidth', 1)
    % Uncomment line below to plot FF
    plot(0:dt:T, rad2deg(FF), 'c.', 'Linewidth', 0.5)
    xlabel('t (s)'), ylabel('$\delta_s$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'delta_', num2str(i),'.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
end

figure(fig_cnt); clf(fig_cnt); hold on;
title(strcat(ES_method, '-ESC | Thrust Coeff'))
plot(0:dt:T, cT, '-', 'Color', [0 0.3 0.9], 'Linewidth', 1.5)
% Uncomment line below to plot filtered cT
% plot(0:dt:T, cT_hat, '--', 'Color', [0.1 0.8 0.8], 'Linewidth', 1.5)
% Uncomment line below to include reference lines
plot(0:dt:T, cT_ref, '--', 'Color', [0.9 0.1 0], 'Linewidth', 1)
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
    xlabel('t (s)'), ylabel('$\hat{\zeta}$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'cT_grad_', num2str(i),'.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
end

if strcmp(ES_method, 'NB')
    figure(fig_cnt); clf(fig_cnt); hold on;
    title('NB-ESC | HPF output')
    plot(0:dt:T, hpf, 'Linewidth', 1.5)
    
    yhpf = zeros(size(hpf));
    for ihpf = 1:N
        yhpf(ihpf) = sum(hpf(1:ihpf)) / ihpf;
    end
    plot(0:dt:T, yhpf, '--')
    
    xlabel('t (s)')
    if save == 1
        filename = fullfile(strcat(dir,'hpf.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1; 

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

    % Hessian reference 
    % Operating point
    hess     = zeros(1, N);
    inv_hess = zeros(1, N);
    dsa      = data.sheeting_angle(2) - data.sheeting_angle(1);

    for k = 1:N
        sa = sheet_angle(:, k+1)';
        
        % Get interpolation axes
        % Axes
        sax = max((sa(1) - dsa), data.sheeting_angle(1)):dsa:min((sa(1) + dsa), data.sheeting_angle(end));

        % Interpolate
        cT_interp = squeeze(interpn(data.AWA, data.sheeting_angle, V, AWA(k), sax));
        
        % Local (numerical) hessian (w/ error condition)        
        if length(cT_interp) < 3
            hess(k)     = nan;
            inv_hess(k) = nan;
        else                
            h = gradient(gradient(cT_interp, sax), sax);
            hess(k) = h(2);
            inv_hess(k) = inv(hess(k));
        end
    end

    figure(fig_cnt); clf(fig_cnt); hold on;
    title('NB-ESC | Hessian Estimate [Averaged]')   
    npoints  = ceil(fs / f_dither);  
    plot(0:dt:T, movmean(squeeze(cT_hessian), npoints))
    % Uncomment line below to plot reference
    plot(0:dt:T, hess, 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}$', 'Interpreter', 'Latex')   
    
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
   
    figure(fig_cnt); clf(fig_cnt); hold on;
    title('NB-ESC | Hessian Inverse Estimate [Averaged]')
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(2:end)), npoints))
    % Uncomment line below to plot reference
    plot(0:dt:T, squeeze(inv_hess(1, 1, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_inv_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    % Inverse of estimated Hessian
    cancel_factor = zeros(size(cT_hessian));
    for k = 1:N
        cancel_factor(k) = cT_hessian_inv(k+1) * cT_hessian(k);
    end

    figure(fig_cnt); clf(fig_cnt); hold on;
    title('NB-ESC | Hessian Cancellation Factor [Averaged]')
    plot(0:dt:T, movmean(squeeze(cancel_factor), npoints))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cancel_factor_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    figure(fig_cnt); clf(fig_cnt); hold on;
    title('NB-ESC | Hessian Cancellation Factor ')
    plot(0:dt:T, squeeze(cancel_factor))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{1, 1}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cancel_factor.fig'));
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

%% Stats
% Error
MSE_sheet_angle   = mean((sheet_angle(:, 2:end) - sheet_angle_ref).^2, 2);
MSE_cT            = mean((cT - cT_ref).^2);

if strcmp(ES_method, 'NB')
    MSE_cancel_factor = mean((movmean(squeeze(cT_hessian_inv(1, 1, 2:end)), npoints)' - 1).^2, 2);
end

% Criterion
cT_accum        = sum(cT, 2);
cT_accum_opt    = sum(cT_ref, 2);

if save == 1
    fileID = fopen(strcat(dir,'diary.txt'),'w');
elseif save == 0
    fileID = 1;
end

% Data string
if strcmp(ES_method, 'GB')
    data_str = sprintf(['Params:\n'...
                        'AWA: %f\n'...
                        'cT Filtering: %s\n'...
                        'cT Filtering param: %f\n'...
                        'HPF: fc = %f\n'...
                        'LPF: fc = %f\n'...
                        'Dithers: A = %s, f = %s\n'...  
                        'Integrator gain (diag): K = %s\n'...
                        'MSE SA: %s\n' ...
                        'MSE cT: %f\n'...
                        'Accumulated cT: %f (optimal = %f)\n'], rad2deg(ship.yaw), cT_filter, cT_filter_param, fc_hp, fc_lp, ...
                                        num2str(A'), num2str(f_dither'), num2str(diag(K)'), num2str(MSE_sheet_angle'), ...
                                        MSE_cT, cT_accum, cT_accum_opt);

elseif strcmp(ES_method, 'NB')
   data_str = sprintf(['Params:\n'...
                        'AWA: %f\n'...
                        'cT Filtering: %s\n'...
                        'cT Filtering param: %f\n'...
                        'HPF: fc = %f\n'...
                        'LPF: fc = %f\n'...
                        'Dithers: A = %s, f = %s\n'...
                        'Ricatti Filter: wric = %f\n'...
                        'Initial Hessian inverse value: %s\n'...
                        'Integrator gain (diag): K = %s\n'...
                        'MSE SA: %s\n' ...
                        'MSE Cancellation Factor : %s\n' ...
                        'MSE cT: %f\n'...
                        'Accumulated cT: %f (optimal = %f)\n'], rad2deg(ship.yaw), cT_filter, cT_filter_param, fc_hp, fc_lp, ...
                                         num2str(A'), num2str(f_dither'), wric, num2str(diag(ric_0)'), ...
                                         num2str(diag(K)'), num2str(MSE_sheet_angle'), num2str(MSE_cancel_factor), ...
                                         MSE_cT, cT_accum, cT_accum_opt);
end

% Print / Save params
fprintf(fileID, "----%s----\n", datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
fprintf(fileID, data_str);
if fileID ~= 1
    fprintf(data_str)
end
if save == 1
    fclose(fileID);
end