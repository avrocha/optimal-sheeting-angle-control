%---
% Gradient- and Newton-based extremum seeking controller for multivariate static maps
% J(\theta) = cT(sheeting_angle)
%---
% Copyright: Alexandre Vieira da Rocha

%% Init
clear; clc; close all;  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

addpath JavaFoil;  addpath Foils;
global ship counter;
fprintf('-------------------------------------------------------------\n');

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

% Criterion
J = @(sheeting_angle)(getfield(calc_objective_mod(sheeting_angle), 'cT'));

% Method
% 'GB' for Gradient based or 'NB' for Newton based
ES_method = 'GB';
% 1 to save figures and diary or 0 to plot figures and print diary
save = 0; 

%% Simulation
fs = 1; % sampling frequency (Hz)
dt = 1/fs;
T  = 200;
N  = length(0:dt:T);
    
ship.yaw = deg2rad(45);
% sheet_angle_0 = [deg2rad(-35); deg2rad(-35)]; % delta(0) [nx1]
sheet_angle_0 = deg2rad(-35); % delta(0) [nx1]

if strcmp(ES_method, 'GB')
    fprintf("Gradient-based ESC selected\n.")

    % Uncomment params below for 1D set of parameters
    f             = 0.1; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc_hp         = 0.09; % HPF cutoff freq
    fc_lp         = 0.09; % LPF cutoff freq
    K             = 0.0750; % gain (>0 since extremum is maximum)
 
    % Uncomment params below for 2D set of parameters
%     f             = [0.15; 0.1]; % dither freq
%     A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
%     fc            = 0.05; % HPF cutoff freq
%     K             = diag([0.0750, 0.0750]); % gain (>0 since extremum is maximum)

    [sheet_angle, cT, cT_grad] = gbesc(J, dt, N, f, A, fc_hp, fc_lp, K, sheet_angle_0, false);

elseif strcmp(ES_method, 'NB') % [WIP]
    fprintf("Newton-based ESC selected\n.")

    % Parameters
%     f             = [0.15; 0.1]; % dither freq
%     A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
%     fc            = 0.05; % HPF cutoff freq
%     K             = diag([0.0025, 0.0025]); % gain (>0 since extremum is maximum)
%     wric          = 2 * pi * 0.05 * min(f); % ricatti filter parameter
%     ric_0         = diag([-30, -30]);
    f             = 0.1; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc_hp         = 0.09; % HPF cutoff freq
    fc_lp         = 0.09; % LPF cutoff freq
    K             = 0.0025; % gain (>0 since extremum is maximum)
    wric          = 2 * pi * 0.03 * f; % ricatti filter parameter
    ric_0         = -30;

    [sheet_angle, cT, cT_grad, cT_hessian, cT_hessian_inv] = nbesc(J, dt, N, f, A, fc_hp, fc_lp, K, sheet_angle_0, wric, ric_0, false);

else
    fprintf("Wrong method. Select either \'GB\' or \'NB\'.\n")
end

%% Plots
fig_cnt = 1;
n = length(sheet_angle_0);

dir = strcat('plots\', ES_method,'_ESC\static_', num2str(n),'D\');

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
                        'Integrator gain (diag): K = %s\n'...
                        'Initial input value (diag): %s\n\n'], rad2deg(ship.yaw), fc_hp, fc_lp, num2str(A'), num2str(f'), num2str(diag(K)'), num2str(sheet_angle_0'));

elseif strcmp(ES_method, 'NB')
   data_str = sprintf(['Params:\n'...
                        'AWA: %f'...
                        'HPF: fc = %f\n'...
                        'LPF: fc = %f\n'...
                        'Dithers: A = %s, f = %s\n'...
                        'Ricatti Filter: wric = %f\n'...
                        'Initial Hessian inverse value: %s\n'...
                        'Integrator gain (diag): K = %s\n'...
                        'Initial input value (diag): %s\n\n'], rad2deg(ship.yaw), fc_hp, fc_lp, num2str(A'), num2str(f'), wric, num2str(diag(ric_0)'), num2str(diag(K)'), num2str(sheet_angle_0));
end

% Print / Save params
fprintf(fileID, "----%s----\n", datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
fprintf(fileID, data_str);
if save == 1
    fclose(fileID);
end
    
% Plots
for i = 1:n
    figure(fig_cnt); clf(fig_cnt); hold on;
    title(strcat(ES_method, '-ESC | Sheeting Angle-', num2str(i)))
    plot(0:dt:T, rad2deg(sheet_angle(i, 1:end-1)), 'b-', 'Linewidth', 2)
    % Uncomment lines below to include reference lines
%     ref = [-20, -30];
%     plot(0:dt:T, -ref(i)*ones(N,1), 'r--', 'Linewidth', 1)
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
end

%% Functions
function [u, y, dy] = gbesc(J, dt, N, f, A, fc_hp, fc_lp, K, u0, lp_bool)
    % Gradient-based extremum seeking controller for 1D static maps
    % Inputs:
    % - J      : optimization criterion [function handle]
    % - dt     : simulation step [s]
    % - N      : simulation lenght 
    % - f      : sinusoidal dithers frequency [Hz]
    % - A      : sinusoidal dithers amplitude [rad]
    % - fc_hp  : HPF cut-off frequency [Hz]
    % - fc_lp  : LPPF cut-off frequency [Hz]
    % - K      : integrator gain
    % - u0     : input initial value
    % - lp_bool: use LPF [boolean] 
    % Outputs:
    % - u : control variable
    % - y : criterion output
    % - dy: criterion gradient estimate
    
    bworder = 5;
    % HPF
    [bh,ah]   = butter(bworder, fc_hp*dt, 'high');
    % LPF
    [bl,al]   = butter(bworder, fc_lp*dt, 'low');

    % Data structures
    n     = length(u0);
    u_hat = [u0, zeros(n, N)];
    u     = [u0, zeros(n, N)];
    y     = zeros(1, N);
    dy    = zeros(n, N);
    hpf   = zeros(1, N);
    zeta  = zeros(1, N);
    lpf   = zeros(1, N);
    
    if (length(A) ~= length(f)) && (length(f) ~= n) && (n ~= length(ship.yaw()))
        fprintf('Dimensions do not match\n.')
        return       
    end

    tic
    for i = 1:N
        if rem(i, 10) == 0
            fprintf("Iteration %d\n", i);
        end
        t = i*dt;
        y(i) = J(u(:, i));
        % Avoid numerical singularities
        if i > 1 && (y(i) > 1.2*y(i-1) || y(i) < 0.8*y(i-1))
                y(i) = y(i-1);
        end

        if i >= bworder+1
            for j = 1:bworder+1
                hpf(i) = hpf(i) + bh(j)*y(i-j+1);
            end
    
            for j = 2:bworder+1
                hpf(i) = hpf(i) - ah(j)*hpf(i-j+1);
            end
            
            hpf(i) = 1/ah(1) * hpf(i);
        else
            for j = 1:i
                hpf(i) = hpf(i) + bh(j)*y(i-j+1);
            end
    
            for j = 2:i
                hpf(i) = hpf(i) - ah(j)*hpf(i-j+1);
            end
            
            hpf(i) = 1/ah(1) * hpf(i);
        end
        
        if lp_bool
            zeta(i) = hpf(i) * sin(2*pi*f*t);
            
            % LPF
            if i >= bworder+1
                for j = 1:bworder+1
                    lpf(i) = lpf(i) + bl(j)*zeta(i-j+1);
                end
        
                for j = 2:bworder+1
                    lpf(i) = lpf(i) - al(j)*lpf(i-j+1);
                end
                
                lpf(i) = 1/al(1) * lpf(i);
            else
                for j = 1:i
                    lpf(i) = lpf(i) + bl(j)*zeta(i-j+1);
                end
        
                for j = 2:i
                    lpf(i) = lpf(i) - al(j)*lpf(i-j+1);
                end
                
                lpf(i) = 1/al(1) * lpf(i);
            end
    
            dy(:, i)      = lpf(i);
        
        else
            dy(:, i) = hpf(i) * sin(2*pi*f*t);

        end

        u_hat(:, i+1) = u_hat(:, i) + dt * K * dy(:, i); % single integrator

        u(:, i+1)     = u_hat(:, i+1) + A .* sin(2*pi*f*t);
        
        % Error condition
        if any(u(i+1) > pi) || any(u(i+1) < -pi) 
            break
        end
    end
    toc
end

function [u, y, dy, ddy, ddy_inv] = nbesc(J, dt, N, f, A, fc_hp, fc_lp, K, u0, wric, ddy0, lp_bool)
    % Gradient-based extremum seeking controller for 1D static maps
    % Inputs:
    % - J      : optimization criterion [function handle]
    % - dt     : simulation step [s]
    % - N      : simulation lenght 
    % - f      : sinusoidal dither frequency [Hz]
    % - A      : sinusoidal dither amplitude [rad]
    % - fc_hp  : HPF cut-off frequency [Hz]
    % - fc_lp  : LPF cut-off frequency [Hz]
    % - K      : integrator gain
    % - u0     : input initial value
    % - wric   : Ricatti filter parameter
    % - ddy0   : hessian inverse initial value
    % - lp_bool: use LPF [boolean] 

    % Outputs:
    % - u      : control variable
    % - y      : criterion output
    % - dy     : criterion gradient estimate
    % - ddy    : hessian estimate
    % - ddy_inv: hessian inverse estimate
    
    % Data structures
    n              = length(u0);
    u_hat          = [u0, zeros(n, N)];
    u              = [u0, zeros(n, N)];
    y              = zeros(1, N);
    dy             = zeros(n, N);
    hpf   = zeros(1, N);
    zeta  = zeros(1, N);
    lpf   = zeros(1, N);
    ddy            = zeros(n, n, N);
    ddy_inv        = zeros(n, n, N+1); % output of ricatti filter
    ddy_inv(:,:,1) = ddy0;

    if (length(A) ~= length(f)) && (length(f) ~= n) && (n ~= length(ship.yaw())) ...
            && (length(ship.yaw()) ~= size(ddy0, 1))
        fprintf('Dimensions do not match\n.')
        return       
    end
    
    bworder = 5;
    % HPF
    [bh,ah]   = butter(bworder, fc_hp*dt, 'high');
    % LPF
    [bl,al]   = butter(bworder, fc_lp*dt, 'low');

    % N(t) - Hessian estimate
    Nhess = cell(n,n);
    for i = 1:n
        for j = 1:n
            if i == j
                Nhess{i,j} = @(t) 16 / (rad2deg(A(i))^2) * (sin(2*pi*f(i)*t)^2 - 0.5);
            else
                Nhess{i,j} = @(t) 4 / (rad2deg(A(i)) * rad2deg(A(j))) * sin(2*pi*f(i)*t) * sin(2*pi*f(j)*t);
            end
        end
    end

    tic
    for i = 1:N
        if rem(i, 10) == 0
            fprintf("Iteration %d\n", i);
        end
        t = i*dt;
        y(i) = J(u(:, i));
        % Avoid numerical singularities
        if i > 1 && (y(i) > 1.2*y(i-1) || y(i) < 0.8*y(i-1))
                y(i) = y(i-1);
        end
        
        % HPF
        if i >= bworder+1
            for j = 1:bworder+1
                hpf(i) = hpf(i) + bh(j)*y(i-j+1);
            end
    
            for j = 2:bworder+1
                hpf(i) = hpf(i) - ah(j)*hpf(i-j+1);
            end
            
            hpf(i) = 1/ah(1) * hpf(i);
        else
            for j = 1:i
                hpf(i) = hpf(i) + bh(j)*y(i-j+1);
            end
    
            for j = 2:i
                hpf(i) = hpf(i) - ah(j)*hpf(i-j+1);
            end
            
            hpf(i) = 1/ah(1) * hpf(i);
        end
        
        if lp_bool
            zeta(i) = hpf(i) * sin(2*pi*f*t);
            
            % LPF
            if i >= bworder+1
                for j = 1:bworder+1
                    lpf(i) = lpf(i) + bl(j)*zeta(i-j+1);
                end
        
                for j = 2:bworder+1
                    lpf(i) = lpf(i) - al(j)*lpf(i-j+1);
                end
                
                lpf(i) = 1/al(1) * lpf(i);
            else
                for j = 1:i
                    lpf(i) = lpf(i) + bl(j)*zeta(i-j+1);
                end
        
                for j = 2:i
                    lpf(i) = lpf(i) - al(j)*lpf(i-j+1);
                end
                
                lpf(i) = 1/al(1) * lpf(i);
            end
    
            dy(:, i)      = lpf(i);
        
        else
            dy(:, i) = hpf(i) * sin(2*pi*f*t);
            
        end

        for j = 1:n
            for k = 1:n
                ddy(j,k,i) = hpf(i) * Nhess{j,k}(t); % \hat{H} = N(t)y_hp
            end
        end

        ddy_inv(:, :, i+1) = ddy_inv(:, :, i) + dt * (wric * ddy_inv(:, :, i) - ...
            wric * ddy_inv(:, :, i) * ddy(:, :, i) * ddy_inv(:, :, i)); % Euler discretization of Ricatti equation    

        u_hat(:, i+1) = u_hat(:, i) - dt * K * ddy_inv(:, :, i+1) * dy(:, i); % single integrator
        
        u(:, i+1) = u_hat(:, i+1) + A .* sin(2*pi*f*t);
        
        if any(u(i+1) > pi) || any(u(i+1) < -pi) 
            break
        end
    end
    toc
end