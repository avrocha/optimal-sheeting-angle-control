%---
% Real-data implementation of ESC controller.
% [Interpolated criterion - Fast computations for 2D scenario]
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
y_abs     = rad2deg(abs(y/n));
y_hat_abs = rad2deg(abs(y_hat/n));

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
R2     = Rig(75,0); % pivot x,y,  
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile, x, y, dx, chord

ship.addRig(R1);
ship.addRig(R2);

ship.yaw = deg2rad(0);
scale    = calc_scale();

% Method
% 'GB' for Gradient based | 'NB' for Newton based
ES_method = 'NB';
% 1 to save figures and diary or 0 to plot figures and print diary
save = 0;

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

% Uncomment lines below for LPF AWA
% Filter AWA - LPF
% [b, a] = butter(5, 0.02/fs, 'low');
% b = fliplr(b);
% a = fliplr(a);
% awa_filt = zeros(1, N);
% for i = 1:N
%     if i > 5
%         awa_filt(i) = 1/a(end) * (-a(1:end-1) * awa_filt(i-6+1:i-1)' + b * AWA(i-6+1:i)');    
%     else 
%         
%         y_init   = [AWA(1) * ones(1, 6-i), AWA(1:i)];
%         awa_filt_init = [AWA(1) * ones(1, 6-i), awa_filt(1:i-1)];            
%         awa_filt(i) = 1/a(end) * (-a(1:end-1) * awa_filt_init' + b * y_init');
%     end
% end
% 
% % Frequency analysis from data
% fx        = (0:N-1)*(fs/N);
% y         = fft(AWA);
% y_filt    = fft(awa_filt);
% y_abs     = rad2deg(abs(y/N));
% y_f_abs   = rad2deg(abs(y_filt/N));
% 
% figure(fig_cnt); clf(fig_cnt); hold on;
% stem(fx, y_abs)
% stem(fx, y_f_abs)
% xlim([0, 0.6])
% ylim([0, 2])
% legend('RAW', 'LPF')
% xlabel('f [Hz]', 'Interpreter', 'Latex')
% ylabel('$|$AWA(jw)$|$ [deg]', 'Interpreter', 'Latex')
% 
% % filename = fullfile('plots\final\AWA_45_freq.eps');
% % saveas(figure(fig_cnt), filename, 'epsc');
% 
% fig_cnt = fig_cnt + 1;
% 
% figure(fig_cnt); clf(fig_cnt); hold on;
% plot(0:dt:T, rad2deg(AWA))
% plot(0:dt:T, rad2deg(awa_filt))
% xlim([0, T])
% legend('RAW', 'LPF')    
% xlabel('t [s]', 'Interpreter', 'latex'), ylabel('AWA [deg]', 'Interpreter', 'latex')
% % 
% % filename = fullfile('plots\final\AWA_45_profile.eps');
% % saveas(figure(fig_cnt), filename, 'epsc');
% 
% fig_cnt = fig_cnt + 1;
%
% AWA = awa_filt; % filtered awa

% % Uncomment lines below for constant AWA
% AWA = deg2rad(100) * ones(1, N);
% AWA = deg2rad(45) * ones(1, N);

% Feedforward (Piecewise constant)
switch data_source
    case 'tacking'
%         Piecewise Linear FF
        FF = zeros(2, N);

        % AWA LPF indexes
%         [~, idx_1] = min(abs(340 - t_sim));
%         [~, idx_2] = min(abs(440 - t_sim));
%         [~, idx_3] = min(abs(1180 - t_sim));
%         [~, idx_4] = min(abs(1240 - t_sim));
%         [~, idx_5] = min(abs(1570 - t_sim));
%         [~, idx_6] = min(abs(1650 - t_sim));

%         AWA raw indexes
        [~, idx_1] = min(abs(311 - t_sim));
        [~, idx_2] = min(abs(350 - t_sim));
        [~, idx_3] = min(abs(1150 - t_sim));
        [~, idx_4] = min(abs(1200 - t_sim));
        [~, idx_5] = min(abs(1540 - t_sim));
        [~, idx_6] = min(abs(1600 - t_sim));
        
        FF(:, 1:idx_1)     = deg2rad(-25);

        m1 = deg2rad(50)/(idx_2-idx_1) * ones(2, 1);
        b1 = deg2rad(25) - m1*idx_2;
        FF(:, idx_1:idx_2) =  m1 .* (idx_1:1:idx_2) + b1;

        FF(:, idx_2:idx_3) = deg2rad(25);
        
        m2 = -deg2rad(50)/(idx_4-idx_3) * ones(2, 1);
        b2 = deg2rad(-25) - m2*idx_4;       
        FF(:, idx_3:idx_4) = m2 .* (idx_3:1:idx_4) + b2;
        
        FF(:, idx_4:idx_5) = deg2rad(-25);
        
        m3 = deg2rad(50)/(idx_6-idx_5) * ones(2, 1);
        b3 = deg2rad(25) - m3*idx_6;
        FF(:, idx_5:idx_6) = m3  .* (idx_5:1:idx_6) + b3;
        
        FF(:, idx_6:end)   = deg2rad(25); 

%         % Piecewise Constant FF
%         FF = zeros(2, N);
%         [~, idx_1] = min(abs(340 - t_sim));
%         [~, idx_2] = min(abs(1180 - t_sim));
%         [~, idx_3] = min(abs(1570 - t_sim));
%         
%         FF(:, 1:idx_1)     = deg2rad(25);
%         FF(:, idx_1:idx_2) = deg2rad(-25);
%         FF(:, idx_2:idx_3) = deg2rad(25);
%         FF(:, idx_3:end) = deg2rad(-25);

        % Constant FF
%         FF = deg2rad(-35) * ones(2, N);
        
    case 'awa_100'
        FF = deg2rad(-85) * ones(2, N);
end

sheet_angle_0 = FF(:, 1);

% Prep Interpolation
% Choose data source
switch data_source
    case 'tacking'
        load('data\measured_data\awa_pm_45\cT_2D.mat')
    case 'awa_100'
        load('data\measured_data\awa_100\cT_2D.mat')
    otherwise
        disp('Error: Select valid data source.\n')
end

X                  = cell(3, 1);
[X{1}, X{2}, X{3}] = ndgrid(data.AWA, data.sheeting_angle_1, data.sheeting_angle_2);
V                  = data.cT;

% Run medfilt2 for every 2D level
for i = 1:size(V, 1)
    Vi        = squeeze(V(i, :, :));
    Vi_median = medfilt2(Vi);
    Vi_diff   = Vi_median - Vi;
    
    Vi(abs(Vi_diff) > 0.1) = Vi_median(abs(Vi_diff) > 0.1);
    V(i, :, :) = Vi;
end

% Filtering cT - {'RAW', 'EMA', 'LPF'}
cT_filter = 'LPF';
switch cT_filter
    case 'EMA' % Param = EMA 'speed'
        switch data_source
            case 'tacking'
                cT_filter_param = 0.2;   
            case 'awa_100'
                cT_filter_param = 0.2;
                
            fprintf("cT filter EMA time constant = %f Hz\n", 1/(log(1-cT_filter_param)*fs))
        end

    case 'LPF' % Param = cut-off frequency for LPF
        cT_filter_param = 4;
        fprintf("cT filter LPF cut off frequency = %f Hz\n", cT_filter_param)

    case 'RAW'
        cT_filter_param = 1;
end

if strcmp(ES_method, 'GB')
    fprintf("Gradient-based ESC selected.\n")
    
    % 2D set of parameters
    f             = 0.01; % tuning param: constant coeff
    delta         = 0.1;  % tuning param: constant coeff

    f_dither      = [20*f; 17*f]; % dither freq (foremost = last array element)
    A             = deg2rad(2) * ones(2, 1); % dither amplitude
    fc_hp         = 3*f; % HPF cutoff freq
    fc_lp         = 1*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    
    % Data-specific tuning
    switch data_source
        case 'tacking'
            ric_0 = eye(2) * -0.5;       
        case 'awa_100'
            ric_0 = eye(2) * -0.05;
    
    end
    
    K = f * delta * 5 * (-ric_0); % gain (>0 since extremum is maximum)

    % Criterion
    Jgb = @(sheeting_angle, ship) (getfield(calc_objective_mod(sheeting_angle, ship), 'cT'));
    % Interpolated criterion
    Jgb_interp = @(sheeting_angle, ship) interp_criterion(X, V, [ship.yaw, sheeting_angle'], 'linear', Jgb, ship);
    
    localShip = ship;
    [GB.sheet_angle, GB.cT, GB.cT_hat, GB.cT_grad] = gbesc(localShip, Jgb_interp, dt, N, f_dither, A, fc_hp, ...
                            fc_lp, K, sheet_angle_0, AWA, lp_bool, cT_filter, cT_filter_param, FF);   

end

if strcmp(ES_method, 'NB')
    fprintf("Newton-based ESC selected.\n")
    
    % 2D set of parameters
    f             = 0.01; % tuning param: constant coeff
    delta         = 0.1;  % tuning param: constant coeff

    f_dither      = [20*f; 17*f]; % dither freq (foremost = last array element)
    A             = deg2rad(2) * ones(2, 1); % dither amplitude
    fc_hp         = 3*f; % HPF cutoff freq
    fc_lp         = 1*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = f * delta * 5 * eye(2); % gain (>0 since extremum is maximum)

    % Data-specific tuning
    switch data_source
        case 'tacking'
            wric  = 2 * pi * (0.005 * f * delta); % ricatti filter parameter
            ric_0 = eye(2) * -0.5;
        case 'awa_100'
            wric  = 2 * pi * (0.01 * f * delta); % ricatti filter parameter 
            ric_0 = eye(2) * -0.05;
    end

    % Criterion 
    Jnb = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 2), 'cT'));    
    % Interpolated criterion - high resolution -> linear interpolation
    Jnb_interp = @(sheeting_angle, ship) interp_criterion(X, V, [ship.yaw, sheeting_angle'], 'linear', Jnb, ship);

    localShip = ship;
    [NB.sheet_angle, NB.cT, NB.cT_hat, NB.cT_grad, NB.cT_hessian, NB.cT_hessian_inv, NB.hpf] = nbesc(localShip, Jnb_interp, dt, N, f_dither, A, ...
        fc_hp, fc_lp, K, sheet_angle_0, wric, ric_0, AWA, lp_bool, cT_filter, cT_filter_param, FF);    
    
end

%% Plots 
close all
switch ES_method
    case 'GB'
        sheet_angle     = GB.sheet_angle;
        cT              = GB.cT;
        cT_hat          = GB.cT_hat;
        cT_grad         = GB.cT_grad;

    case 'NB'
        sheet_angle     = NB.sheet_angle;
        cT              = NB.cT;
        cT_hat          = NB.cT_hat;
        cT_grad         = NB.cT_grad;
        cT_hessian      = NB.cT_hessian;
        cT_hessian_inv  = NB.cT_hessian_inv;
        hpf             = NB.hpf;
end

n = length(f_dither);

% Directory to save pictures in
dir = ['plots\final\AWA_100\2D\', ES_method,'_ESC\'];

% Check directory
if ~exist(dir, 'dir')
    mkdir(dir)
end

% References
% Get reference in original grid
sa_ref = zeros(2, length(data.AWA));
cT_ref = zeros(1, length(data.AWA));

for i = 1:length(data.AWA)
    [cT_ref(i), I] = max(squeeze(V(i, :, :)), [], 'all');    
    n_max = find(squeeze(V(i, :, :)) == cT_ref(i));
    [row, col]     = ind2sub(size(V, [2, 3]), I);
    sa_ref(1, i)   = data.sheeting_angle_1(row);
    sa_ref(2, i)   = data.sheeting_angle_2(col);

%     fprintf("REF GENERATION - DEBUG\n")
%     fprintf("AWA = %f | cT_opt = %f | SA_opt(1) = %f | SA_opt(2) = %f | n_max = %d\n", ...
%                 rad2deg(data.AWA(i)), cT_ref(i), rad2deg(sa_ref(1, i)), ...
%                 rad2deg(sa_ref(2, i)), length(n_max));
end

% Interpolate references
sheet_angle_ref       = zeros(2, length(AWA));
sheet_angle_ref(1, :) = interp1(data.AWA, sa_ref(1, :), AWA);
sheet_angle_ref(2, :) = interp1(data.AWA, sa_ref(2, :), AWA);
cT_ref                = interp1(data.AWA, cT_ref, AWA);

% Plots
for i = 1:n
    figure(fig_cnt); clf(fig_cnt); hold on;
    title(strcat(ES_method, '-ESC | Sheeting Angle-', num2str(i)))
    plot(0:dt:T, rad2deg(sheet_angle(i, 1:end-1)), 'b-', 'Linewidth', 2)
%     Uncomment line below to include reference lines
    plot(0:dt:T, rad2deg(sheet_angle_ref(i, :)), 'r--', 'Linewidth', 1)
    % Uncomment line below to plot FF
%     plot(0:dt:T, rad2deg(FF), 'c.', 'Linewidth', 0.5)
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
plot(0:dt:T, cT_hat, '--', 'Color', [0.1 0.8 0.8], 'Linewidth', 1.5)
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
    hess     = zeros(2, 2, N);
    inv_hess = zeros(2, 2, N);
    dsa1     = data.sheeting_angle_1(2) - data.sheeting_angle_1(1);
    dsa2     = data.sheeting_angle_2(2) - data.sheeting_angle_2(1);

    for k = 1:N
        sa = sheet_angle(:, k+1)';
        
        % Get interpolation axes
        % Axes
        sa1 = max((sa(1) - dsa1), data.sheeting_angle_1(1)):dsa1:min((sa(1) + dsa1), data.sheeting_angle_1(end));
        sa2 = max((sa(2) - dsa2), data.sheeting_angle_2(1)):dsa2:min((sa(2) + dsa2), data.sheeting_angle_2(end));

        % Interpolate
        cT_interp = squeeze(interpn(data.AWA, data.sheeting_angle_1, data.sheeting_angle_2, V, ...
                                AWA(k), sa1, sa2));
        
        % Local (numerical) hessian (w/ error condition)        
        if any(size(cT_interp) < 3)
            hess(:, :, k)     = NaN(2);
            inv_hess(:, :, k) = NaN(2);
        else                
            [gx, gy] = gradient(cT_interp, sa1, sa1); 
            [gxx, gxy] = gradient(gx, sa1, sa1);
            [~, gyy] = gradient(gy, sa1, sa1);
            
            % Results
            hess(:, :, k) = [gxx(2, 2), gxy(2, 2);
                             gxy(2, 2), gyy(2, 2)];
            
            inv_hess(:, :, k) = inv(hess(:, :, k));
        end
    end

    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Estimate [Averaged]')
    
    coeffs_f = int16(f_dither ./ f);    
    F_common = gcd(coeffs_f(1), coeffs_f(2));
    T_common = 1 / (double(F_common)*f);
    npoints  = T_common * fs;
    
    subplot(3, 1, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(1,1,:)), npoints))
    % Uncomment line below to plot reference
    plot(0:dt:T, squeeze(hess(1,1,:)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{1, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 1, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(2, 2,:)), npoints))
    % Uncomment line below to plot reference
    plot(0:dt:T, squeeze(hess(2, 2,:)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{2, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 1, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(2, 1,:)), npoints))
    % Uncomment line below to plot reference
    plot(0:dt:T, squeeze(hess(2, 1,:)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{2, 1} = \hat{H}_{1, 2}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
   
    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Inverse Estimate [Averaged]')
    
    subplot(3, 1, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(1,1, 2:end)), npoints))
    % Uncomment line below to plot reference
    plot(0:dt:T, squeeze(inv_hess(1, 1, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{1, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 1, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(2, 2, 2:end)), npoints))
    % Uncomment line below to plot reference
    plot(0:dt:T, squeeze(inv_hess(2, 2, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{2, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 1, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(2, 1, 2:end)), npoints))
    % Uncomment line below to plot reference
    plot(0:dt:T, squeeze(inv_hess(2, 1,:)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{2, 1} = \Gamma_{1, 2}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_inv_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    % Inverse of estimated Hessian
    cancel_factor = zeros(size(cT_hessian));
    for k = 1:N
        cancel_factor(:, :, k) = cT_hessian_inv(:, :, k+1) * cT_hessian(:, :, k);
    end

    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Cancellation Factor [Averaged]')
    
    subplot(3, 1, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(1, 1, :)), npoints))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{1, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 1, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(2, 2, :)), npoints))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{2, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 1, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(2, 1, :)), npoints))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{2, 1} = (\Gamma H)_{1, 2}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cancel_factor_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Cancellation Factor ')
    
    subplot(3, 1, 1); hold on;
    plot(0:dt:T, squeeze(cancel_factor(1, 1, :)))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{1, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 1, 2); hold on;
    plot(0:dt:T, squeeze(cancel_factor(2, 2, :)))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{2, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 1, 3); hold on;
    plot(0:dt:T, squeeze(cancel_factor(2, 1, :)))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{2, 1} = (\Gamma H)_{1, 2}$', 'Interpreter', 'Latex')
    
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
%     saveas(figure(fig_cnt), filename);
end

%% Stats
% Error
MSE_sheet_angle   = mean((sheet_angle(:, 2:end) - sheet_angle_ref).^2, 2);
MSE_cT            = mean((cT - cT_ref).^2);

if strcmp(ES_method, 'NB')
    MSE_cancel_factor    = zeros(1, 3);
    MSE_cancel_factor(1) = mean((movmean(squeeze(cT_hessian_inv(1,1, 2:end)), npoints)' - 1).^2, 2);
    MSE_cancel_factor(2) = mean((movmean(squeeze(cT_hessian_inv(2,2, 2:end)), npoints)' - 1).^2, 2);
    MSE_cancel_factor(3) = mean((movmean(squeeze(cT_hessian_inv(2,1, 2:end)), npoints)' - 1).^2, 2);
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


%% Comparison plots - Used in report/presentation
%---------------------------------
% Plots for final presentation/report
% - Run main section once for GB and once for NB, without clearing the
% workspace
%---------------------------------
if ~exist('GB', 'var') || ~exist('NB', 'var')
    disp('Error: Missing either GB or NB results.')
    return
end
close all
save = 0;

switch data_source
    case 'tacking'
        dir = 'plots\final\AWA_45\2D\report\';

    case 'awa_100'
        dir = 'plots\final\AWA_100\2D\report\';
end

% Check directory
if ~exist(dir, 'dir')
    mkdir(dir)
end

n = length(f_dither);

% References
% Get reference in original grid
sa_ref = zeros(2, length(data.AWA));
cT_ref = zeros(1, length(data.AWA));

for i = 1:length(data.AWA)
    [cT_ref(i), I] = max(squeeze(V(i, :, :)), [], 'all');    
    n_max = find(squeeze(V(i, :, :)) == cT_ref(i));
    [row, col]     = ind2sub(size(V, [2, 3]), I);
    sa_ref(1, i)   = data.sheeting_angle_1(row);
    sa_ref(2, i)   = data.sheeting_angle_2(col);
end

% Interpolate references
sheet_angle_ref       = zeros(2, length(AWA));
sheet_angle_ref(1, :) = interp1(data.AWA, sa_ref(1, :), AWA);
sheet_angle_ref(2, :) = interp1(data.AWA, sa_ref(2, :), AWA);
cT_ref                = interp1(data.AWA, cT_ref, AWA);

if AWA(1) == deg2rad(100) && ~any(diff(AWA) == 1) 
    sheet_angle_ref(1, :) = deg2rad(-77.4322);
    sheet_angle_ref(2, :) = deg2rad(-80.3182);  
    cT_ref                = 1.151 * ones(1, N);
elseif AWA(1) == deg2rad(45) && ~any(diff(AWA) == 1) 
    sheet_angle_ref(1, :) = deg2rad(-19.637);
    sheet_angle_ref(2, :) = deg2rad(-31.453);  
    cT_ref                = 0.826267 * ones(1, N);
end

% Plots

% SHEETING ANGLES 
% ---
for i = 1:n
    figure(fig_cnt); clf(fig_cnt); hold on;
    plot(0:dt:T, rad2deg(GB.sheet_angle(i, 1:end-1)), '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.5)
    plot(0:dt:T, rad2deg(NB.sheet_angle(i, 1:end-1)), '-', 'color', [0.8500 0.3250 0.0980], 'Linewidth', 1.5)
%     Uncomment line below to include reference lines
    plot(0:dt:T, rad2deg(sheet_angle_ref(i, :)), '--', 'color', [0.8100 0.8100 0.2700],  'Linewidth', 0.5)
    % Uncomment line below to plot FF
    plot(0:dt:T, rad2deg(FF(i, :)), 'c--', 'Linewidth', 1.2)
    xlabel('t [s]', 'Interpreter', 'Latex')
    if i == 1
        ylabel('$\delta_1$ [deg]', 'Interpreter', 'Latex')
        legend('GBESC', 'NBESC', '$\delta_1^*$', 'FF', 'Interpreter', 'Latex', 'location', 'best')
    else
        ylabel('$\delta_2$ [deg]', 'Interpreter', 'Latex')
        legend('GBESC', 'NBESC', '$\delta_2^*$', 'FF', 'Interpreter', 'Latex', 'location', 'best')    
    end

    if save == 1
        filename = fullfile(strcat(dir,'delta_', num2str(i),'.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
end

% CRITERION cT
% ---
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, GB.cT, '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.3)
plot(0:dt:T, NB.cT, '-', 'color', [0.8500 0.3250 0.0980], 'Linewidth', 1.3)
% Uncomment line below to include filtered cT
% plot(0:dt:T, NB.cT_hat, '--', 'Color', [0.1 0.8 0.8], 'Linewidth', 1.5)
% Uncomment line below to include reference lines
plot(0:dt:T, cT_ref, '--', 'color', [0.8100 0.8100 0.2700], 'Linewidth', 0.5)
xlabel('t [s]', 'Interpreter', 'Latex'), ylabel('$J$', 'Interpreter', 'Latex')
legend('GBESC', 'NBESC', '$J^*$', 'Interpreter', 'Latex', 'location', 'best')
if save == 1
    filename = fullfile(strcat(dir,'cT.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% Common period
coeffs_f = int16(f_dither ./ f);    
F_common = gcd(coeffs_f(1), coeffs_f(2));
T_common = 1 / (double(F_common)*f);
npoints  = T_common * fs;

% CANCEL FACTOR H*GAMMA
% ---
cancel_factor = zeros(size(NB.cT_hessian));
for k = 1:N
    cancel_factor(:, :, k) = NB.cT_hessian_inv(:, :, k+1) * NB.cT_hessian(:, :, k);
end

figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(cancel_factor(1, 1, :)), npoints), '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(2, 1, :)), npoints), '-', 'color', [0.8100 0.8100 0.2700], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(2, 2, :)), npoints), '-', 'color', [0.9290 0.6940 0.1250], 'Linewidth', 1.3)
plot(0:dt:T, ones(1, N), 'k--')
plot(0:dt:T, zeros(1, N), 'k--')
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$(\hat{H}\Gamma)_{1,1}$', '$(\hat{H}\Gamma)_{1, 2} = (\hat{H}\Gamma)_{2, 1}$', '$(\hat{H}\Gamma)_{2, 2}$', '','','Interpreter', 'Latex', 'location', 'best')
xlabel('t (s)')
if save == 1
    filename = fullfile(strcat(dir,'cancel_factor_avg.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% HESSIAN ESTIMATE (NBESC)
% ---
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(1, 1, :)), npoints), '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(1, 2, :)), npoints), '-', 'color', [0.8100 0.8100 0.2700], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(2, 2, :)), npoints), '-', 'color', [0.9290 0.6940 0.1250], 'Linewidth', 1.3)
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$\hat{H}_{1,1}$', '$\hat{H}_{1, 2} = \hat{H}_{2, 1}$', '$\hat{H}_{2, 2}$', 'Interpreter', 'Latex', 'location', 'best')
if save == 1
    filename = fullfile(strcat(dir,'cT_hessian.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% HESSIAN INVERSE ESTIMATE (NBESC)
% ---
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(1, 1, 2:end)), npoints), '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(1, 2, 2:end)), npoints), '-', 'color', [0.8100 0.8100 0.2700], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(2, 2, 2:end)), npoints), '-', 'color', [0.9290 0.6940 0.1250], 'Linewidth', 1.3)
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$\Gamma_{1,1}$', '$\Gamma_{1, 2} = \Gamma_{2, 1}$', '$\Gamma_{2, 2}$', 'Interpreter', 'Latex', 'location', 'best')
if save == 1
    filename = fullfile(strcat(dir,'cT_hessian_inv.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% STATS
% --- 
GB.MAE_sheet_angle   = mean(abs(GB.sheet_angle(:, 2:end) - sheet_angle_ref), 2);
NB.MAE_sheet_angle   = mean(abs(NB.sheet_angle(:, 2:end) - sheet_angle_ref), 2);
GB.MAE_cT            = mean(abs(GB.cT - cT_ref));
NB.MAE_cT            = mean(abs(NB.cT - cT_ref));

% Criterion
GB.cT_accum  = sum(GB.cT, 2);
NB.cT_accum  = sum(NB.cT, 2);
cT_accum_opt = sum(cT_ref, 2);

if save == 1
    fileID = fopen(strcat(dir,'diary.txt'),'w');
elseif save == 0
    fileID = 1;
end

data_str = sprintf(['Params'...
                    '\tHPF: fc = %f\n'...
                    '\tLPF (%d): fc = %f\n'...
                    '\tDithers: A = %s, f = %s\n'...
                    '\tRicatti Filter: wric = %f\n'...
                    '\tInitial Hessian inverse value: %s\n'...
                    '\tIntegrator gain (diag): K = %s\n'...
                    'GB\n'...
                    '\tMAE SA: %s\n' ...
                    '\tMAE cT: %f\n'...
                    '\tAccumulated cT: %f\n'...
                    'NB\n'...
                    '\tMAE SA: %s\n' ...
                    '\tMAE cT: %f\n'...
                    '\tAccumulated cT: %f\n'...
                    'Optimum\n'...
                    'Accumulated cT: %f\n'], fc_hp, lp_bool, fc_lp, num2str(A'), num2str(f_dither'), wric, num2str(diag(ric_0)'), ...
                                             num2str(diag(K)'),num2str(GB.MAE_sheet_angle'), GB.MAE_cT, GB.cT_accum,...
                                             num2str(NB.MAE_sheet_angle'), NB.MAE_cT, NB.cT_accum,...
                                             cT_accum_opt);

% Print / Save params
fprintf(fileID, "----%s----\n", datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
fprintf(fileID, data_str);
if fileID ~= 1
    fprintf(data_str)
end

if save == 1
    fclose(fileID);
end