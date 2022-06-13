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

% Filter data - LPF, EMA, ORIG
filter_type = 'LPF';

switch filter_type 
    case 'LPF'
        disp('AWA: Butterworth LPF filter.')
        % IIR design
        bworder  = 5;
        fc       = 0.05;
        dt       = 1 / fs_data;
        [b, a]   = butter(bworder, fc*dt, 'low');
        M = bworder + 1; % filter length
        
        % Flip coeffs for vector-wise multiplication
        a = fliplr(a);
        b = fliplr(b);
        
        % Add prefix
        x = [awa(1)*ones(1, bworder+1), awa'];
        y = [awa(1)*ones(1, bworder+1), zeros(1, n)];
        
        for i = M:n+M
            y(i) = 1/a(end) * (-a(1:end-1) * y(i-M+1:i-1)' + b * x(i-M+1:i)');
        end
        
        % Delete initial prefix
        y = y(M+1:end);

    case 'EMA'
        disp('AWA: Exponential moving average (EMA) filter.')
        % Exponential moving average
        y      = [awa(1), zeros(1, n-1)];
        alpha  = 0.005; 
        
        for i = 2:n
            y(i) = (1-alpha)*y(i-1) + alpha*awa(i);
        end

    case 'ORIG'
        disp('AWA: Original filter.')
        y = awa_hat;

    case 'RAW'
        disp('AWA: Raw data.')

    otherwise
        disp('Wrong AWA filter selection.')

end

if ~strcmp(filter_type, 'RAW')
    figure(fig_cnt); clf(fig_cnt);
    hold on;
    plot(1:n, rad2deg(awa'));
    plot(1:n, rad2deg(y));
    ylabel('AWA [rad]')
    xlabel('t [s]')
    legend('RAW', filter_type);
    fig_cnt = fig_cnt + 1;
    
    % Output final AWA to next section
    awa = y;
else
    figure(fig_cnt); clf(fig_cnt);
    hold on;
    plot(1:n, rad2deg(awa'));
    ylabel('AWA [rad]')
    xlabel('t [s]')
    fig_cnt = fig_cnt + 1;
end

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
ES_method = 'NB';
% 1 to save figures and diary or 0 to plot figures and print diary
save = 1;

% Simulation
fs = 10; 
dt = 1/fs;
T  = ceil(time(end));
N  = length(0:dt:T);

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
        sheet_angle_0 = deg2rad(-85);
end

% Prep Interpolation
% Choose data source
load('data\measured_data\awa_100\cT_1D.mat')
%     load('data\measured_data\awa_pm_45\cT_1D.mat')
X            = cell(2, 1);
[X{1}, X{2}] = ndgrid(data.AWA, data.sheeting_angle);
V            = data.cT;

if strcmp(ES_method, 'GB')
    fprintf("Gradient-based ESC selected.\n")
    
    % 1D set of parameters
    f             = 0.01; % tuning param: constant coeff
    delta         = 0.1;  % tuning param: constant coeff

    f_dither      = 20*f; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc_hp         = 2*f; % HPF cutoff freq
    fc_lp         = 15*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = -f * delta * 0.7 * (-45); % gain (>0 since extremum is maximum)
 
    % Criterion
    Jgb = @(sheeting_angle, ship) (getfield(calc_objective_mod(sheeting_angle, ship), 'cT'));
    % Interpolated criterion
    Jgb_interp = @(sheeting_angle, ship) interp_criterion(X, V, [ship.yaw, sheeting_angle'], 'linear', Jgb, ship);
    
    localShip = ship;
    [sheet_angle, cT, cT_grad] = gbesc(localShip, Jgb_interp, dt, N, f_dither, A, fc_hp, fc_lp, K, sheet_angle_0, AWA, lp_bool);   

end

if strcmp(ES_method, 'NB')
    fprintf("Newton-based ESC selected.\n")

    % 1D set of parameters
    f             = 0.01; % tuning param: constant coeff
    delta         = 0.1;  % tuning param: constant coeff

    f_dither      = 20*f; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc_hp         = 2*f; % HPF cutoff freq
    fc_lp         = 15*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = f * delta * 0.7; % gain (>0 since extremum is maximum)
    wric          = 2 * pi * (1 * f * delta); % ricatti filter parameter
    ric_0         = -45;
    
    % Criterion
    Jnb = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 2), 'cT'));    
    % Interpolated criterion - high resolution -> linear interpolation
    Jnb_interp = @(sheeting_angle, ship) interp_criterion(X, V, [ship.yaw, sheeting_angle'], 'linear', Jnb, ship);

    localShip = ship;
    [sheet_angle, cT, cT_grad, cT_hessian, cT_hessian_inv] = nbesc(localShip, Jnb_interp, dt, N, f_dither, A, fc_hp, fc_lp, K, sheet_angle_0, wric, ric_0, AWA, lp_bool);    
    
end

%% Plots 
dir = ['plots\7m_data_AWA_100\real\', ES_method,'_ESC\', filter_type,'_awa\f_0_20\'];

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
plot(0:dt:T,  cT_ref, 'r--', 'Linewidth', 1)
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

%% Stats
MSE_sheet_angle = mean((sheet_angle(2:end) - sheet_angle_ref).^2);
MSE_cT          = mean((cT - cT_ref).^2);
cT_accum        = sum(cT, 2);

if save == 1
    fileID = fopen(strcat(dir,'diary.txt'),'w');
elseif save == 0
    fileID = 1;
end


% Data string
if strcmp(ES_method, 'GB')
    data_str = sprintf(['Params:\n'...
                        'AWA: %f\n'...
                        'AWA Filtering: %s\n'...
                        'HPF: fc = %f\n'...
                        'LPF: fc = %f\n'...
                        'Dithers: A = %s, f = %s\n'...  
                        'Integrator gain (diag): K = %s\n'...
                        'MSE SA: %f\n' ...
                        'MSE cT: %f\n'...
                        'Accumulated cT: %f\n'], rad2deg(ship.yaw), filter_type, fc_hp, fc_lp, num2str(A'), ...
                                         num2str(f_dither'), num2str(diag(K)'), MSE_sheet_angle, ...
                                         MSE_cT, cT_accum);

elseif strcmp(ES_method, 'NB')
   data_str = sprintf(['Params:\n'...
                        'AWA: %f\n'...
                        'AWA Filtering: %s\n'...
                        'HPF: fc = %f\n'...
                        'LPF: fc = %f\n'...
                        'Dithers: A = %s, f = %s\n'...
                        'Ricatti Filter: wric = %f\n'...
                        'Initial Hessian inverse value: %s\n'...
                        'Integrator gain (diag): K = %s\n'...
                        'MSE SA: %f\n' ...
                        'MSE cT: %f\n'...
                        'Accumulated cT: %f\n'], rad2deg(ship.yaw), filter_type, fc_hp, fc_lp, num2str(A'), ...
                                         num2str(f_dither'), wric, num2str(diag(ric_0)'), ...
                                         num2str(diag(K)'), MSE_sheet_angle, MSE_cT, cT_accum);
end

% Print / Save params
fprintf(fileID, "----%s----\n", datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
fprintf(fileID, data_str);
fprintf(data_str)
if save == 1
    fclose(fileID);
end