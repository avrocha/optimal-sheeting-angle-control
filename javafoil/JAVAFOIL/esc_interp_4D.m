%---
% Real-data implementation of ESC controller.
% [Interpolated criterion - Fast computations for 4D scenario]
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

% Menus.
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
% 
% [~, itime_max] = min(abs(310 - time));
% time    = time(1:itime_max);
% awa     = out(1:itime_max, 2);
% awa_hat = out(1:itime_max, 3);
% heading = out(1:itime_max, 4);

% Small deviation from 5Hz
fs_data = 1/(time(2) - time(1));
n       = size(out, 1);
% n       = size(time, 1);

% Edit time to obtain constant sampling frequency from the sensor
fs_data = 5; % Both datasets are approximately 5 Hz
time    = (0:1/fs_data:(1/fs_data)*(n-1))';

% Frequency analysis from data
fx        = (0:n-1)*(fs_data/n);
y         = fft(awa);
y_hat     = fft(awa_hat);
y_abs     = abs(y/n);
y_hat_abs = abs(y_hat/n);

% Uncomment lines below to plot frequency content of AWA data
figure(fig_cnt); clf(fig_cnt);
subplot(2, 1, 1); hold on;
plot(fx, rad2deg(y_abs));
title('Frequency spectrum of AWA raw data')
xlabel('f [Hz]')
ylabel('|awa(jw)|')

subplot(2, 1, 2); hold on;
plot(fx, rad2deg(y_hat_abs));
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
R3     = Rig(122,0); % pivot x,y,  
R3.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
R4     = Rig(170,0); % pivot x,y,  
R4.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

ship.addRig(R1);
ship.addRig(R2);
ship.addRig(R3);
ship.addRig(R4);

ship.yaw = deg2rad(0);
scale    = calc_scale();

% Method
% 'GB' for Gradient based | 'NB' for Newton based
ES_method = 'GB';
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
[b, a] = butter(5, 0.02/fs, 'low');
b = fliplr(b);
a = fliplr(a);
awa_filt = zeros(1, N);
for i = 1:N
    if i > 5
        awa_filt(i) = 1/a(end) * (-a(1:end-1) * awa_filt(i-6+1:i-1)' + b * AWA(i-6+1:i)');    
    else         
        y_init   = [AWA(1) * ones(1, 6-i), AWA(1:i)];
        awa_filt_init = [AWA(1) * ones(1, 6-i), awa_filt(1:i-1)];            
        awa_filt(i) = 1/a(end) * (-a(1:end-1) * awa_filt_init' + b * y_init');
    end
end

% Frequency analysis from data
fx        = (0:N-1)*(fs/N);
y         = fft(AWA);
y_filt    = fft(awa_filt);
y_abs     = rad2deg(abs(y/N));
y_f_abs   = rad2deg(abs(y_filt/N));

figure(fig_cnt); clf(fig_cnt); hold on;
stem(fx(1:3:N), y_abs(1:3:N))
stem(fx(1:3:N), y_f_abs(1:3:N))
xlim([0, 0.6])
ylim([0, 5])
legend('RAW', 'LPF')
xlabel('f [Hz]', 'Interpreter', 'Latex')
ylabel('$|$AWA(jw)$|$ [deg]', 'Interpreter', 'Latex')
% 
% filename = fullfile('plots\final\AWA_45_freq_both.eps');
% saveas(figure(fig_cnt), filename, 'epsc');

fig_cnt = fig_cnt + 1;

figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, rad2deg(AWA))
% plot(0:dt:T, rad2deg(awa_filt))
xlim([0, T])
legend('RAW', 'LPF')    
xlabel('t [s]', 'Interpreter', 'latex'), ylabel('AWA [deg]', 'Interpreter', 'latex')

filename = fullfile('plots\final\AWA_100.eps');
% saveas(figure(fig_cnt), filename, 'epsc');

fig_cnt = fig_cnt + 1;

AWA      = awa_filt; % filtered awa

% % Uncomment lines below for constant AWA
% AWA = deg2rad(100) * ones(1, N);
% AWA = deg2rad(45) * ones(1, N);

% Feedforward (Piecewise constant)
switch data_source
    case 'tacking'
        % Piecewise Linear FF
%         FF = zeros(4, N);
%         [~, idx_1] = min(abs(310 - t_sim));
%         [~, idx_2] = min(abs(1150 - t_sim));
%         [~, idx_3] = min(abs(1560 - t_sim));
%         
%         FF(:, 1:idx_1) = [deg2rad(-20);
%                           deg2rad(-30);
%                           deg2rad(-35);
%                           deg2rad(-40)];
% 
%         FF(:, idx_1:idx_2) = [deg2rad(30);
%                               deg2rad(35);
%                               deg2rad(38);
%                               deg2rad(43)];
% 
%         FF(:, idx_2:idx_3) = [deg2rad(-20);
%                               deg2rad(-30);
%                               deg2rad(-35);
%                               deg2rad(-40)];
% 
%         FF(:, idx_3:end) = [deg2rad(30);
%                               deg2rad(35);
%                               deg2rad(38);
%                               deg2rad(43)]; 
        
%         FF = [deg2rad(-20) * ones(1, N);
%               deg2rad(-30) * ones(1, N);
%               deg2rad(-35) * ones(1, N);
%               deg2rad(-40) * ones(1, N)]; 
        
        % Constant FF
        FF = [deg2rad(-5) * ones(1, N);
              deg2rad(-20) * ones(1, N);
              deg2rad(-25) * ones(1, N);
              deg2rad(-30) * ones(1, N)]; 

    case 'awa_100'
        FF = [deg2rad(-82) * ones(1, N);
              deg2rad(-84) * ones(1, N);
              deg2rad(-86) * ones(1, N);
              deg2rad(-87) * ones(1, N)];

end

sheet_angle_0 = FF(:, 1);

% Prep Interpolation
% Choose data source
switch data_source
    case 'tacking'
        load('data\measured_data\awa_pm_45\cT_4D_test.mat')
    case 'awa_100'
        load('data\measured_data\awa_100\cT_4D.mat')
    otherwise
        disp('Error: Select valid data source.\n')
end

V = data.cT;

for i = 1:length(data.AWA)
    % Run median filter on data
    Vi      = squeeze(V(i, :, :, :, :));
    Vi_med  = medfilt4(squeeze(V(i, :, :, :, :)));
    Vi_diff = abs(Vi_med - Vi);
    
    Vi(Vi_diff > 0.05) = Vi_med(Vi_diff > 0.05);
    
    V(i, :, :, :, :) = Vi;
end

% Filtering cT - {'RAW', 'EMA', 'LPF'}
cT_filter = 'RAW';
switch cT_filter
    case 'EMA' % Param = EMA 'speed'
        switch data_source
            case 'tacking'
                cT_filter_param = 0.1;   
            case 'awa_100'
                cT_filter_param = 0.1;
                
            fprintf("cT filter EMA time constant = %f Hz\n", 1/(log(1-cT_filter_param)*fs))
        end

    case 'LPF' % Param = cut-off frequency for LPF
        cT_filter_param = 1;
        fprintf("cT filter LPF cut off frequency = %f Hz\n", cT_filter_param)

    case 'RAW'
        cT_filter_param = 1;
end

if strcmp(ES_method, 'GB')
    fprintf("Gradient-based ESC selected.\n")
    
    % 4D set of parameters
    f             = 0.01; % tuning param: constant coeff
    delta         = 0.1;  % tuning param: constant coeff

    f_dither      = [14*f; 16*f; 18*f; 20*f]; % dither freq (foremost = last array element)
    A             = deg2rad(1)*ones(4, 1); % dither amplitude
    fc_hp         = 3*f; % HPF cutoff freq
    fc_lp         = 5*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    
    switch data_source
        case 'tacking'
            ric_0 = diag(diag([-0.3921    0.0066    0.0189    0.0443;
                                 0.0066   -0.2149   -0.0018    0.0013;
                                 0.0189   -0.0018   -0.1591   -0.0052;
                                 0.0443    0.0013   -0.0052   -0.0932])); % Obtained with simulated AWA
            
            ric_0 = eye(4) * -0.2; % For simulated AWA
        
        case 'awa_100'            
           ric_0 = diag(diag([ -0.4972   -0.0389   -0.0152    0.0514;
                               -0.0389   -0.2209   -0.0129   -0.0033;
                               -0.0152   -0.0129   -0.1744   -0.0022;
                               0.0514   -0.0033   -0.0022   -0.1013])); % Obtained with simulated AWA
            
           ric_0 = eye(4) * -0.2; % For simulated AWA
    end
    
    K = f * delta * 300 * (-ric_0); % gain (>0 since extremum is maximum)    
    
    % Criterion
    Jgb = @(sheeting_angle, ship) (getfield(calc_objective_mod(sheeting_angle, ship), 'cT'));
    % Interpolated criterion
    Jgb_interp = @(sheeting_angle, ship) interp_criterion_irregular(data.AWA, data.sheeting_angle, V, [ship.yaw, sheeting_angle'], 'linear', Jgb, ship);
    localShip = ship;
    [GB.sheet_angle, GB.cT, GB.cT_hat, GB.cT_grad, GB.hpf] = gbesc(localShip, Jgb_interp, dt, N, f_dither, A, fc_hp, ...
                            fc_lp, K, sheet_angle_0, AWA, lp_bool, cT_filter, cT_filter_param, FF);   

end

if strcmp(ES_method, 'NB')
    fprintf("Newton-based ESC selected.\n")

    % 4D set of parameters
    f             = 0.01; % tuning param: constant coeff
    delta         = 0.1;  % tuning param: constant coeff
    f_dither      = [8*f; 14*f; 18*f; 20*f]; % dither freq (foremost = last array element)
%     f_dither      = [0.4; 0.7; 0.9; 1];
    A             = deg2rad(1) * ones(4, 1); % dither amplitude
    fc_hp         = 3*f; % HPF cutoff freq
    fc_lp         = 5*f; % LPF cutoff freq
    lp_bool       = false; % Use LPF
    K             = f * delta * 150 * eye(4); % gain (>0 since extremum is maximum)

    switch data_source
        case 'tacking'
            wric  = 2 * pi * (0.1 * f * delta); % ricatti filter parameter

            ric_0 = diag(diag([-0.3921    0.0066    0.0189    0.0443;
                                 0.0066   -0.2149   -0.0018    0.0013;
                                 0.0189   -0.0018   -0.1591   -0.0052;
                                 0.0443    0.0013   -0.0052   -0.0932])); % Obtained with simulated AWA
            
            ric_0 = eye(4) * -0.2; % For simulated AWA
        
        case 'awa_100'
            wric  = 2 * pi * (1.5 * f * delta); % ricatti filter parameter 
            
            ric_0 = diag(diag([ -0.4972   -0.0389   -0.0152    0.0514;
                               -0.0389   -0.2209   -0.0129   -0.0033;
                               -0.0152   -0.0129   -0.1744   -0.0022;
                               0.0514   -0.0033   -0.0022   -0.1013])); % Obtained with simulated AWA
            
            ric_0 = eye(4) * -0.2; % For simulated AWA
    end
        
    % Criterion 
    Jnb = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 2), 'cT'));    
    % Interpolated criterion - high resolution -> linear interpolation
    Jnb_interp = @(sheeting_angle, ship) interp_criterion_irregular(data.AWA, data.sheeting_angle, V, [ship.yaw, sheeting_angle'], 'linear', Jnb, ship);
    
    localShip = ship;
    [NB.sheet_angle, NB.cT, NB.cT_hat, NB.cT_grad, NB.cT_hessian, NB.cT_hessian_inv, NB.hpf] = nbesc(localShip, Jnb_interp, dt, N, f_dither, A, ...
        fc_hp, fc_lp, K, sheet_angle_0, wric, ric_0, AWA, lp_bool, cT_filter, cT_filter_param, FF);    
    
end

%% Plots 
switch ES_method
    case 'GB'
        sheet_angle     = GB.sheet_angle;
        cT              = GB.cT;
        cT_hat          = GB.cT_hat;
        cT_grad         = GB.cT_grad;
        hpf             = GB.hpf;

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

dir = ['plots\final\AWA_45_sim\4D\', ES_method,'_ESC\'];

% Check directory
if ~exist(dir, 'dir')
    mkdir(dir)
end

n = length(f_dither);

% References
% Get reference in original grid
sa_ref = zeros(4, length(data.AWA));
cT_ref = zeros(1, length(data.AWA));

for i = 1:length(data.AWA)
    [cT_ref(i), I] = max(squeeze(V(i, :, :, :, :)), [], 'all', 'linear');    
    idxs           = find(squeeze(V(i, :, :, :, :)) == cT_ref(i));

    [idx1, idx2, idx3, idx4] = ind2sub(size(V, 2:5), I);
    sa_ref(1, i)             = data.sheeting_angle(1, idx1, i);
    sa_ref(2, i)             = data.sheeting_angle(2, idx2, i);
    sa_ref(3, i)             = data.sheeting_angle(3, idx3, i);
    sa_ref(4, i)             = data.sheeting_angle(4, idx4, i);
end 

% Interpolate references
sheet_angle_ref       = zeros(4, length(AWA));
sheet_angle_ref(1, :) = interp1(data.AWA, sa_ref(1, :), AWA);
sheet_angle_ref(2, :) = interp1(data.AWA, sa_ref(2, :), AWA);
sheet_angle_ref(3, :) = interp1(data.AWA, sa_ref(3, :), AWA);
sheet_angle_ref(4, :) = interp1(data.AWA, sa_ref(4, :), AWA);
cT_ref                = interp1(data.AWA, cT_ref, AWA);

% Plots
for i = 1:n
    figure(fig_cnt); clf(fig_cnt); hold on;
    title(strcat(ES_method, '-ESC | Sheeting Angle-', num2str(i)))
    plot(0:dt:T, rad2deg(sheet_angle(i, 1:end-1)), 'b-', 'Linewidth', 2)
    % Uncomment line below to include reference lines
    plot(0:dt:T, rad2deg(sheet_angle_ref(i, :)), 'r--', 'Linewidth', 1)
    % Uncomment line below to plot FF
    plot(0:dt:T, rad2deg(FF(i, :)), 'c', 'Linewidth', 0.5)
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

figure(fig_cnt); clf(fig_cnt); hold on;
title(strcat(ES_method, '-ESC | HPF'))
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
    xlabel('t (s)'), ylabel('$\Gamma$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_inv.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    % Hessian reference - Low numerical accuracy
    hess     = zeros(4, 4, N);
    inv_hess = zeros(4, 4, N);
%     
%     for k = 1:N
%         % Get AWA neighbors
%         [~, idx_neigh] = min(abs(data.AWA - AWA(k)));
%         if data.AWA(idx_neigh) > AWA(k)
%             neighbors = [idx_neigh - 1, idx_neigh];
%         elseif data.AWA(idx_neigh) < AWA(k)
%             neighbors = [idx_neigh, idx_neigh + 1];
%         else
%             neighbors = idx_neigh;
%         end
%         
%         % Interpolate cT neighborhood in SA
%         cT_interp = zeros(length(neighbors), 3, 3, 3, 3);
%         sa_axes   = zeros(4, 3);
%         dsa       = data.sheeting_angle(1, 2, 1) - data.sheeting_angle(1, 1, 1); % Common resolution to the 4 axes
%         i = 1;
%         for idx = neighbors
%             sa_axes(1, :) = max(sheet_angle_0(1)-dsa, data.sheeting_angle(1, 1, idx)):dsa:min(sheet_angle_0(1)+dsa, data.sheeting_angle(1, end, idx));
%             sa_axes(2, :) = max(sheet_angle_0(2)-dsa, data.sheeting_angle(2, 1, idx)):dsa:min(sheet_angle_0(2)+dsa, data.sheeting_angle(2, end, idx));
%             sa_axes(3, :) = max(sheet_angle_0(3)-dsa, data.sheeting_angle(3, 1, idx)):dsa:min(sheet_angle_0(3)+dsa, data.sheeting_angle(3, end, idx));
%             sa_axes(4, :) = max(sheet_angle_0(4)-dsa, data.sheeting_angle(4, 1, idx)):dsa:min(sheet_angle_0(4)+dsa, data.sheeting_angle(4, end, idx));
%             
%             [X1, X2, X3, X4]     = ndgrid(squeeze(data.sheeting_angle(1, :, idx)), squeeze(data.sheeting_angle(2, :, idx)), ...
%                                     squeeze(data.sheeting_angle(3, :, idx)), squeeze(data.sheeting_angle(4, :, idx)));
%             [Xq1, Xq2, Xq3, Xq4] = ndgrid(sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
%         
%             interp_res = squeeze(interpn(X1, X2, X3, X4, squeeze(data.cT(idx, :, :, :, :)), ...
%                                     Xq1, Xq2, Xq3, Xq4)); 
%                 
%             if any(size(interp_res) < 3)
%                 disp('Hess Ref: Boundary problems in interpolation.')
%                 cT_interp(i, :, :, :, :) = NaN(3, 3, 3, 3);
%             else
%                 cT_interp(i, :, :, :, :) = interp_res;
%             end
%         
%             i = i + 1;
%         
%         end
%         
%         % Interpolate neighborhood in AWA
%         cT_interp2 = zeros(3, 3, 3, 3);
%         if size(cT_interp, 1) > 1
%             cT_1 = squeeze(cT_interp(1, :, :, :, :));
%             cT_2 = squeeze(cT_interp(2, :, :, :, :));
%             for i = 1:3^4
%                 cT_interp2(i) = interp1(data.AWA(neighbors), [cT_1(i), cT_2(i)], AWA(k));
%             end
%         else
%             cT_interp2 = cT_interp;
%         end
%         
%         % Get Hessian
%         [g1, g2, g3, g4]     = gradient(cT_interp2, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
%         [g11, g12, g13, g14] = gradient(g1, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
%         [g21, g22, g23, g24] = gradient(g2, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
%         [g31, g32, g33, g34] = gradient(g3, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
%         [g41, g42, g43, g44] = gradient(g4, sa_axes(1, :), sa_axes(2, :), sa_axes(3, :), sa_axes(4, :));
%         
%         hess(:, :, k)     = [g11(2), g12(2), g13(2), g14(2);
%                                 g21(2), g22(2), g23(2), g24(2);
%                                 g31(2), g32(2), g33(2), g34(2);
%                                 g41(2), g42(2), g43(2), g44(2)];
%                     
%         inv_hess(:, :, k) = inv(hess(:, :, k));
%     end

    coeffs_f = int16(f_dither ./ f);    
    F_common = gcd(coeffs_f(1), gcd(coeffs_f(2), gcd(coeffs_f(3), coeffs_f(4))));
    T_common = 1 / (double(F_common)*f);
    npoints  = T_common * fs;
    
    % Hessian Estimate
    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Estimate [Averaged] (Diag)')
    
    subplot(4, 1, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(1,1,:)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(1,1,:)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{1, 1}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(2, 2,:)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(2, 2,:)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{2, 2}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(3, 3,:)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(3, 3, :)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{3, 3}$', 'Interpreter', 'Latex')

    subplot(4, 1, 4); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(4, 4,:)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(4, 4, :)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{4, 4}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_diag_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Estimate [Averaged] (Cross)')
    
    subplot(3, 2, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(2, 1,:)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(2, 1, :)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{2, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(3, 1,:)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(3, 1,:)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{3, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(3, 2, :)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(3, 2, :)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{3, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 4); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(4, 1, :)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(4, 1, :)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{4, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 5); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(4, 2, :)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(4, 2, :)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{4, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 6); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian(4, 3, :)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(hess(4, 3, :)), 'r--')
    xlabel('t (s)'), ylabel('$\hat{H}_{4, 3}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_cross_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
   
    % Hessian Inverse estimate
    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Inverse Estimate [Averaged] (Diag)')
    
    subplot(4, 1, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(1, 1, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(1,1,:)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{1, 1}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(2, 2, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(2, 2,:)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{2, 2}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(3, 3, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(3, 3, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{3, 3}$', 'Interpreter', 'Latex')

    subplot(4, 1, 4); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(4, 4, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(4, 4, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{4, 4}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_inv_diag_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Inverse Estimate [Averaged] (Cross)')
    
    subplot(3, 2, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(2, 1, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(2, 1, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{2, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(3, 1, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(3, 1,:)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{3, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(3, 2, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(3, 2, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{3, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 4); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(4, 1, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(4, 1, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{4, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 5); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(4, 2, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(4, 2, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{4, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 6); hold on;
    plot(0:dt:T, movmean(squeeze(cT_hessian_inv(4, 3, 2:end)), npoints))
    % Uncomment line below to plot reference
%     plot(0:dt:T, squeeze(inv_hess(4, 3, :)), 'r--')
    xlabel('t (s)'), ylabel('$\Gamma_{4, 3}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_inv_cross_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;
    
    % Cancellation factor 
    cancel_factor = zeros(size(cT_hessian));
    for k = 1:N
        cancel_factor(:, :, k) = cT_hessian_inv(:, :, k+1) * cT_hessian(:, :, k);
    end

    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Cancellation Factor [Averaged] (Diag)')
    
    subplot(4, 1, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(1, 1, :)), npoints))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{1, 1}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(2, 2, :)), npoints))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{2, 2}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(3, 3, :)), npoints))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{3, 3}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 4); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(3, 3, :)), npoints))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{4, 4}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cancel_factor_diag_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Cancellation Factor [Averaged] (Cross)')
    
    subplot(3, 2, 1); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(2, 1, :)), npoints))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{2, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 2); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(3, 1, :)), npoints))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{3, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 3); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(3, 2, :)), npoints))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{3, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 4); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(4, 1, :)), npoints))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{4, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 5); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(4, 2, :)), npoints))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{4, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 6); hold on;
    plot(0:dt:T, movmean(squeeze(cancel_factor(4, 3, :)), npoints))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{4, 3}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cancel_factor_cross_avg.fig'));
        saveas(figure(fig_cnt), filename);
    end
   
    fig_cnt = fig_cnt + 1;

    % Cancellation factor - Not averaged
    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Cancellation Factor (Diag)')
    
    subplot(4, 1, 1); hold on;
    plot(0:dt:T, squeeze(cancel_factor(1, 1, :)))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{1, 1}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 2); hold on;
    plot(0:dt:T, squeeze(cancel_factor(2, 2, :)))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{2, 2}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 3); hold on;
    plot(0:dt:T, squeeze(cancel_factor(3, 3, :)))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{3, 3}$', 'Interpreter', 'Latex')
    
    subplot(4, 1, 4); hold on;
    plot(0:dt:T, squeeze(cancel_factor(3, 3, :)))
    plot(0:dt:T, ones(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{4, 4}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cancel_factor_diag.fig'));
        saveas(figure(fig_cnt), filename);
    end
    fig_cnt = fig_cnt + 1;

    figure(fig_cnt); clf(fig_cnt); hold on;
    sgtitle('NB-ESC | Hessian Cancellation Factor (Cross)')
    
    subplot(3, 2, 1); hold on;
    plot(0:dt:T, squeeze(cancel_factor(2, 1, :)))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{2, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 2); hold on;
    plot(0:dt:T, squeeze(cancel_factor(3, 1, :)))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{3, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 3); hold on;
    plot(0:dt:T, squeeze(cancel_factor(3, 2, :)))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{3, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 4); hold on;
    plot(0:dt:T, squeeze(cancel_factor(4, 1, :)))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{4, 1}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 5); hold on;
    plot(0:dt:T, squeeze(cancel_factor(4, 2, :)))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{4, 2}$', 'Interpreter', 'Latex')
    
    subplot(3, 2, 6); hold on;
    plot(0:dt:T, squeeze(cancel_factor(4, 3, :)))
    plot(0:dt:T, zeros(1, N), 'r--', 'Linewidth', 0.8)
    xlabel('t (s)'), ylabel('$(\Gamma H)_{4, 3}$', 'Interpreter', 'Latex')
    
    if save == 1
        filename = fullfile(strcat(dir,'cancel_factor_cross.fig'));
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
MSE_sheet_angle = mean((sheet_angle(:, 2:end) - sheet_angle_ref).^2, 2);
MSE_cT          = mean((cT - cT_ref).^2);
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
                        'MSE cT: %f\n'...
                        'Accumulated cT: %f (optimal = %f)\n'], rad2deg(ship.yaw), cT_filter, cT_filter_param, fc_hp, fc_lp, ...
                                         num2str(A'), num2str(f_dither'), wric, num2str(diag(ric_0)'), ...
                                         num2str(diag(K)'), num2str(MSE_sheet_angle'), MSE_cT, cT_accum, cT_accum_opt);
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
%---------------------------------
if ~exist('GB', 'var') || ~exist('NB', 'var')
    disp('Error: Missing either GB or NB results.')
    return
end
close all
save = 0;

switch data_source
    case 'tacking'
%         dir = 'plots\final\AWA_45_sim\4D\report\';
        dir = 'plots\final\AWA_45_filt\4D\report\';

    case 'awa_100'
%         dir = 'plots\final\AWA_100_sim\4D\report\';
        dir = 'plots\final\AWA_100_filt\4D\report\higher_gb_f\';
%         dir = 'plots\final\AWA_100\4D\report\';
end

% Check directory
if ~exist(dir, 'dir')
    mkdir(dir)
end

n = length(f_dither);

% References
% Get reference in original grid
sa_ref = zeros(4, length(data.AWA));
cT_ref = zeros(1, length(data.AWA));

for i = 1:length(data.AWA)
    [cT_ref(i), I] = max(squeeze(V(i, :, :, :, :)), [], 'all');    
    n_max          = find(squeeze(V(i, :, :, :, :)) == cT_ref(i));

    [idx1, idx2, idx3, idx4] = ind2sub(size(V, 2:5), I);
    sa_ref(1, i)             = data.sheeting_angle(1, idx1, i);
    sa_ref(2, i)             = data.sheeting_angle(2, idx2, i);
    sa_ref(3, i)             = data.sheeting_angle(3, idx3, i);
    sa_ref(4, i)             = data.sheeting_angle(4, idx4, i);

    fprintf("REF GENERATION - DEBUG\n")
    fprintf("AWA = %f | cT_opt = %f | SA_opt(1) = %f | SA_opt(2) = %f | SA_opt(3) = %f | SA_opt(4) = %f | n_max = %d\n", ...
                rad2deg(data.AWA(i)), cT_ref(i), rad2deg(sa_ref(1, i)), ...
                rad2deg(sa_ref(2, i)), rad2deg(sa_ref(3, i)), ...
                rad2deg(sa_ref(4, i)), length(n_max));
end 

% Interpolate references
sheet_angle_ref       = zeros(4, length(AWA));
sheet_angle_ref(1, :) = interp1(data.AWA, sa_ref(1, :), AWA);
sheet_angle_ref(2, :) = interp1(data.AWA, sa_ref(2, :), AWA);
sheet_angle_ref(3, :) = interp1(data.AWA, sa_ref(3, :), AWA);
sheet_angle_ref(4, :) = interp1(data.AWA, sa_ref(4, :), AWA);
cT_ref                = interp1(data.AWA, cT_ref, AWA);

if AWA(1) == deg2rad(100) && ~any(diff(AWA) == 1) 
    sheet_angle_ref(1, :) = deg2rad(-74.8244);
    sheet_angle_ref(2, :) = deg2rad(-77.1475);  
    sheet_angle_ref(3, :) = deg2rad(-79.1362);  
    sheet_angle_ref(4, :) = deg2rad(-81.1574);  
    cT_ref                = 0.978 * ones(1, N);
elseif AWA(1) == deg2rad(45) && ~any(diff(AWA) == 1) 
    sheet_angle_ref(1, :) = deg2rad(-13.7);
    sheet_angle_ref(2, :) = deg2rad(-23.822);  
    sheet_angle_ref(3, :) = deg2rad(-29.9);  
    sheet_angle_ref(4, :) = deg2rad(-35.2);  
    cT_ref                = 0.743 * ones(1, N);
end

% Plots

% SHEETING ANGLES 
% ---
for i = 1:n
    figure(fig_cnt); clf(fig_cnt); hold on;
    plot(0:dt:T, rad2deg(GB.sheet_angle(i, 1:end-1)), '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.5)
    plot(0:dt:T, rad2deg(NB.sheet_angle(i, 1:end-1)), '-', 'color', [0.8500 0.3250 0.0980], 'Linewidth', 1.5)
%     Uncomment line below to include reference lines
    plot(0:dt:T, rad2deg(sheet_angle_ref(i, :)), 'k--', 'Linewidth', 0.5)
    % Uncomment line below to plot FF
%     plot(0:dt:T, rad2deg(FF(i, :)), 'c--', 'Linewidth', 1.2)
    xlabel('t [s]', 'Interpreter', 'Latex')
    if i == 1
        ylabel('$\delta_1$ [deg]', 'Interpreter', 'Latex')
        legend('GBESC', 'NBESC', '$\delta_1^*$', 'FF', 'Interpreter', 'Latex', 'location', 'best')
    elseif i == 2
        ylabel('$\delta_2$ [deg]', 'Interpreter', 'Latex')
        legend('GBESC', 'NBESC', '$\delta_2^*$', 'FF', 'Interpreter', 'Latex', 'location', 'best')    
    elseif i == 3
        ylabel('$\delta_3$ [deg]', 'Interpreter', 'Latex')
        legend('GBESC', 'NBESC', '$\delta_3^*$', 'FF', 'Interpreter', 'Latex', 'location', 'best')    
    elseif i == 4
        ylabel('$\delta_4$ [deg]', 'Interpreter', 'Latex')
        legend('GBESC', 'NBESC', '$\delta_4^*$', 'FF', 'Interpreter', 'Latex', 'location', 'best')    
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
plot(0:dt:T, cT_ref, 'k--', 'Linewidth', 0.5)
xlabel('t [s]', 'Interpreter', 'Latex'), ylabel('$J$', 'Interpreter', 'Latex')
legend('GBESC', 'NBESC', '$J^*$', 'Interpreter', 'Latex', 'location', 'best')
if save == 1
    filename = fullfile(strcat(dir,'cT.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% Common period
coeffs_f = int16(f_dither ./ f);    
F_common = gcd(coeffs_f(1), gcd(coeffs_f(2), gcd(coeffs_f(3), coeffs_f(4))));
T_common = 1 / (double(F_common)*f);
npoints  = T_common * fs;

% CANCEL FACTOR H*GAMMA
% ---
cancel_factor = zeros(size(NB.cT_hessian));
for k = 1:N
    cancel_factor(:, :, k) = NB.cT_hessian_inv(:, :, k+1) * NB.cT_hessian(:, :, k);
end

% Diagonal terms
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(cancel_factor(1, 1, :)), npoints), '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(2, 2, :)), npoints), '-', 'color', [0.8100 0.8100 0.2700], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(3, 3, :)), npoints), '-', 'color', [0.9290 0.6940 0.1250], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(4, 4, :)), npoints), '-', 'color', [0.6350 0.0780 0.1840], 'Linewidth', 1.3)
plot(0:dt:T, ones(1, N), 'k--')
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$(\hat{H}\Gamma)_{1,1}$', '$(\hat{H}\Gamma)_{2,2}$', '$(\hat{H}\Gamma)_{3,3}$', '$(\hat{H}\Gamma)_{4,4}$', '', 'Interpreter', 'Latex', 'location', 'best')
xlabel('t (s)')
if save == 1
    filename = fullfile(strcat(dir,'cancel_factor_diag.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% Cross terms
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(cancel_factor(2, 1, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(3, 1, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(3, 2, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(4, 1, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(4, 2, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(cancel_factor(4, 3, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, zeros(1, N), 'k--')
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$(\hat{H}\Gamma)_{2, 1}$', '$(\hat{H}\Gamma)_{3, 1}$', '$(\hat{H}\Gamma)_{3, 2}$', '$(\hat{H}\Gamma)_{3, 2}$', '$(\hat{H}\Gamma)_{4, 1}$', ...
    '$(\hat{H}\Gamma)_{4, 2}$', '$(\hat{H}\Gamma)_{4, 3}$', '', 'Interpreter', 'Latex', 'location', 'best')
xlabel('t (s)')
if save == 1
    filename = fullfile(strcat(dir,'cancel_factor_cross.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% HESSIAN ESTIMATE (NBESC)
% ---
% Diagonal terms
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(1, 1, :)), npoints), '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(2, 2, :)), npoints), '-', 'color', [0.8100 0.8100 0.2700], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(3, 3, :)), npoints), '-', 'color', [0.9290 0.6940 0.1250], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(4, 4, :)), npoints), '-', 'color', [0.6350 0.0780 0.1840], 'Linewidth', 1.3)
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$\hat{H}_{1,1}$', '$\hat{H}_{2,2}$', '$\hat{H}_{3,3}$', '$\hat{H}_{4,4}$', '','Interpreter', 'Latex', 'location', 'best')
xlabel('t (s)')
if save == 1
    filename = fullfile(strcat(dir,'cT_hessian_diag.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% Cross terms
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(2, 1, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(3, 1, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(3, 2, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(4, 1, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(4, 2, :)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian(4, 3, :)), npoints), '-', 'Linewidth', 1.3)
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$\hat{H}_{2, 1}$', '$\hat{H}_{3, 1}$', '$\hat{H}_{3, 2}$', '$\hat{H}_{3, 2}$', '$\hat{H}_{4, 1}$', ...
    '$\hat{H}_{4, 2}$', '$\hat{H}_{4, 3}$', 'Interpreter', 'Latex')
xlabel('t (s)')
if save == 1
    filename = fullfile(strcat(dir,'cT_hessian_cross.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% HESSIAN INVERSE ESTIMATE (NBESC)
% ---
% Diagonal terms
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(1, 1, 2:end)), npoints), '-', 'color', [0 0.4470 0.7410], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(2, 2, 2:end)), npoints), '-', 'color', [0.8100 0.8100 0.2700], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(3, 3, 2:end)), npoints), '-', 'color', [0.9290 0.6940 0.1250], 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(4, 4, 2:end)), npoints), '-', 'color', [0.6350 0.0780 0.1840], 'Linewidth', 1.3)
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$\Gamma_{1,1}$', '$\Gamma_{2,2}$', '$\Gamma_{3,3}$', '$\Gamma_{4,4}$', '','Interpreter', 'Latex', 'location', 'best')
xlabel('t (s)')
if save == 1
    filename = fullfile(strcat(dir,'cT_hessian_inv_diag.fig'));
    saveas(figure(fig_cnt), filename);
end
fig_cnt = fig_cnt + 1;

% Cross terms
figure(fig_cnt); clf(fig_cnt); hold on;
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(2, 1, 2:end)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(3, 1, 2:end)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(3, 2, 2:end)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(4, 1, 2:end)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(4, 2, 2:end)), npoints), '-', 'Linewidth', 1.3)
plot(0:dt:T, movmean(squeeze(NB.cT_hessian_inv(4, 3, 2:end)), npoints), '-', 'Linewidth', 1.3)
xlabel('t [s]', 'Interpreter', 'Latex')
legend('$\Gamma_{2, 1}$', '$\Gamma_{3, 1}$', '$\Gamma_{3, 2}$', '$\Gamma_{3, 2}$', '$\Gamma_{4, 1}$', ...
    '$\Gamma_{4, 2}$', '$\Gamma_{4, 3}$', 'Interpreter', 'Latex')
xlabel('t (s)')
if save == 1
    filename = fullfile(strcat(dir,'cT_hessian_inv_cross.fig'));
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