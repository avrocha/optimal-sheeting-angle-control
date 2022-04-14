%---
% Gradient- and Newton-based extremum seeking controller for a 1D static map
% Future work: adapt functions for {2,4}D static maps
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
ship.addRig(R1);
ship.yaw = deg2rad(0);
scale    = calc_scale();

% Criterion
J = @(sheeting_angle)(getfield(calc_objective_mod(sheeting_angle), 'cT'));

% Method
% 'GB' for Gradient based or 'NB' for Newton based
ES_method = 'NB';
% 1 to save figures and diary or 0 to plot figures and print diary
save = 0; 

%% Simulation
fs = 1; % sampling frequency (Hz)
dt = 1/fs;
T  = 1000;
N = length(0:dt:T);

AWA = deg2rad(45) + deg2rad(30)*sin(2*pi/T * (0:dt:T)) + deg2rad(normrnd(0,1,1,N)); 
sheet_angle_0 = deg2rad(0);

if strcmp(ES_method, 'GB')
    fprintf("Gradient-based ESC selected\n.")

    % Parameters
    f             = 0.1; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc            = 0.05; % HPF cutoff freq
    K             = 0.075; % gain (>0 since extremum is maximum)

    [sheet_angle, cT, cT_grad] = gbesc_1d(J, dt, N, f, A, fc, K, sheet_angle_0, AWA);

elseif strcmp(ES_method, 'NB')
    fprintf("Newton-based ESC selected\n.")

    % Parameters
    f             = 0.1; % dither freq
    A             = deg2rad(2); % dither amplitude
    fc            = 0.05; % HPF cutoff freq
    K             = 0.0025; % gain (>0 since extremum is maximum)
    Ahess         = 16 / (rad2deg(A)^2);  % hessian estimate amplitude
    fhess         = f; % hessian estimate frequency
    wric          = 2 * pi * 0.05 * f; % ricatti filter parameter
    ric_0         = -30;
    
    [sheet_angle, cT, cT_grad, cT_hessian, cT_hessian_inv] = nbesc_1d(J, dt, N, f, A, fc, K, ...
                                    sheet_angle_0, Ahess, fhess, wric, ric_0, AWA);

else
    fprintf("Wrong method. Select either \'GB\' or \'NB\'.\n")
    
end

%% Plots
dir = strcat('plots\',ES_method,'_ESC\time_variant\dist\');
if save == 1
    fileID = fopen(strcat(dir,'diary.txt'),'a');
elseif save == 0
    fileID = 1;
end

% Get reference from structure
    load('cT_SA_AWA.mat')
    sheet_angle_ref = zeros(length(AWA), 1);
    cT_ref          = zeros(length(AWA), 1);
    [~, iy_max] = min(abs(data.y_grid(1, :) - 20));
    for i = 1:length(AWA)
        [~, ix] = min(abs(data.x_grid(:, 1) - rad2deg(AWA(i))));
        [cT_ref(i), iy] = max(data.cT(ix, 1:iy_max));
        sheet_angle_ref(i) = data.y_grid(1,iy);
    end

if strcmp(ES_method, 'GB')
    % Print/Save params
    fprintf(fileID, "----%s----\n", datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
    fprintf(fileID,['Params:\n'...
                    'AWA limits: [%f, %f]\n'...
                    'HPF: fc = %f\n'...
                    'Main dither: A = %f, f = %d\n'...  
                    'Integrator gain: K = %f\n'...
                    'Initial input value: %f\n\n'], [min(rad2deg(AWA)); max(rad2deg(AWA)); fc; A; f; K; sheet_angle_0]);
    if save == 1
        fclose(fileID);
    end
    
    figure(1); clf(1); hold on;
    title('GB-ESC | Sheeting Angle')
    plot(0:dt:T, rad2deg(sheet_angle(1:end-1)), 'b-', 'Linewidth', 2)
    plot(0:dt:T, sheet_angle_ref, 'r--', 'Linewidth', 1.3)
    xlabel('t [s]'), ylabel('deg [º]')
    legend('$\delta_s$', '$\delta_s^*$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'delta.fig'));
        saveas(figure(1), filename);
    end

    figure(2); clf(2); hold on;
    title('GB-ESC | Thrust Coeff')
    plot(0:dt:T, cT, 'b-', 'Linewidth', 2)
    plot(0:dt:T, cT_ref, 'r--', 'Linewidth', 1.3)
    xlabel('t [s]')
    legend('cT', 'cT opt')
    if save == 1
        filename = fullfile(strcat(dir,'cT.fig'));
        saveas(figure(2), filename);
    end

    figure(3); clf(3); hold on;
    title('GB-ESC | Gradient Estimate $\hat{\zeta}$', 'Interpreter', 'Latex')
    plot(0:dt:T, cT_grad, 'b-', 'Linewidth', 2)
    plot(0:dt:T, movmean(cT_grad, 10), 'r--', 'Linewidth', 1.3)
    xlabel('t [s]')
    if save == 1
        filename = fullfile(strcat(dir,'cT_grad.fig'));
        saveas(figure(3), filename);
    end

    figure(4); clf(4); 
    title('GB-ESC | AWA')
    plot(0:dt:T, rad2deg(AWA), 'b-', 'Linewidth', 2)
    xlabel('t (s)'), ylabel('deg (º)')
    if save == 1
        filename = fullfile(strcat(dir,'AWA.fig'));
        saveas(figure(4), filename);
    end

elseif strcmp(ES_method, 'NB')
    % Print/Save params
    fprintf(fileID, "----%s----\n", datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
    fprintf(fileID,['Params:\n'...
                    'AWA limits: [%f, %f]\n'...
                    'HPF: fc = %f\n'...
                    'Main dither: A = %f, f = %d\n'...
                    'Ricatti Filter: wric = %f\n'...
                    'Hessian estimate dither: Ahess =%f, whess = %f\n'...
                    'Integrator gain: K = %f\n'...
                    'Initial input value: %f\n\n'], [min(rad2deg(AWA)); max(rad2deg(AWA)); fc; A; f; wric; Ahess; fhess; K; sheet_angle_0]);
    if save == 1
        fclose(fileID);
    end
    
    figure(1); clf(1); hold on;
    title('NB-ESC | Sheeting Angle')
    plot(0:dt:T, rad2deg(sheet_angle(1:end-1)), 'b-', 'Linewidth', 2)
    plot(0:dt:T, sheet_angle_ref, 'r--', 'Linewidth', 1.3)
    xlabel('t [s]'), ylabel('deg [º]')
    legend('$\delta_s$', '$\delta_s^*$', 'Interpreter', 'Latex')
    if save == 1
        filename = fullfile(strcat(dir,'delta.fig'));
        saveas(figure(1), filename);
    end
    
    figure(2); clf(2); hold on;
    title('NB-ESC | Thrust Coeff')
    plot(0:dt:T, cT, 'b-', 'Linewidth', 2)
    plot(0:dt:T, cT_ref, 'r--', 'Linewidth', 1.3)
    xlabel('t [s]')
    legend('cT', 'cT opt')
    if save == 1
        filename = fullfile(strcat(dir,'cT.fig'));
        saveas(figure(2), filename);
    end
    
    figure(3); clf(3); hold on;
    title('NB-ESC | Gradient Estimate $\hat{\zeta}$', 'Interpreter', 'Latex')
    plot(0:dt:T, cT_grad, 'b-', 'Linewidth', 2)
    plot(0:dt:T, movmean(cT_grad, 10), 'r--', 'Linewidth', 1.3)
    xlabel('t [s]')
    if save == 1
        filename = fullfile(strcat(dir,'cT_grad.fig'));
        saveas(figure(3), filename);
    end
    
    figure(4); clf(4); hold on;
    title('NB-ESC | Hessian Estimate $\hat{H}$', 'Interpreter', 'Latex')
    plot(0:dt:T, cT_hessian, 'b-', 'Linewidth', 2)
    plot(0:dt:T, movmean(cT_hessian, 10), 'r--', 'Linewidth', 1.3)
    xlabel('t [s]')
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian.fig'));
        saveas(figure(4), filename);
    end

    figure(5); clf(5); hold on;
    title('NB-ESC | Hessian Inverse Estimate $\hat{H}^{-1}$', 'Interpreter', 'Latex')
    plot(0:dt:T, cT_hessian_inv(2:end), 'b-', 'Linewidth', 2)
    plot(0:dt:T, movmean(cT_hessian_inv(2:end), 10), 'r--', 'Linewidth', 1.3)
    xlabel('t [s]')
    if save == 1
        filename = fullfile(strcat(dir,'cT_hessian_inv.fig'));
        saveas(figure(5), filename);
    end
    
    figure(6); clf(6); 
    title('NB-ESC | AWA')
    plot(0:dt:T, rad2deg(AWA), 'b-', 'Linewidth', 2)
    xlabel('t [s]'), ylabel('deg [º]')
    if save == 1
        filename = fullfile(strcat(dir,'AWA.fig'));
        saveas(figure(6), filename);
    end
end

%% Functions
function [u, y, dy] = gbesc_1d(J, dt, N, f, A, fc, K, u0, AWA)
    % Gradient-based extremum seeking controller for 1D static maps
    % Inputs:
    % - J  : optimization criterion [function handle]
    % - dt : simulation step [s]
    % - N  : simulation lenght 
    % - f  : sinusoidal dither frequency [Hz]
    % - A  : sinusoidal dither amplitude [rad]
    % - fc : HPF cut-off frequency [Hz]
    % - K  : integrator gain
    % - u0 : input initial value
    % - AWA: time variant AWA
    % Outputs:
    % - u : control variable
    % - y : criterion output
    % - dy: criterion gradient estimate
   
    global ship;

    % HPF
    bworder = 2;
    [b,a]   = butter(bworder, fc*dt, 'high');

    % Data structures
    u_hat = [u0; zeros(N, 1)];
    u     = [u0; zeros(N, 1)];
    y     = zeros(N, 1);
    dy    = zeros(N, 1);
    hpf   = zeros(N, 1);
    
    tic
    for i = 1:N
        if rem(i, 10) == 0
            fprintf("Iteration %d\n", i);
        end
        t = i*dt;
   
        ship.yaw = AWA(i);
        y(i) = J(u(i));
        
        if i >= bworder+1
            for j = 1:bworder+1
                hpf(i) = hpf(i) + b(j)*y(i-j+1);
            end
    
            for j = 2:bworder+1
                hpf(i) = hpf(i) - a(j)*hpf(i-j+1);
            end
            
            hpf(i) = 1/a(1) * hpf(i);
        end
        
        dy(i) = hpf(i) * sin(2*pi*f*t);
        
        u_hat(i+1) = u_hat(i) + K * dt * dy(i); % single integrator
            
        u(i+1) = u_hat(i+1) + A * sin(2*pi*f*t);
        
        % Error condition
        if u(i+1) > pi || u(i+1) < -pi 
            break
        end
    end
    toc
end

function [u, y, dy, ddy, ddy_inv] = nbesc_1d(J, dt, N, f, A, fc, K, u0, Ahess, fhess, wric, ddy0, AWA)
    % Gradient-based extremum seeking controller for 1D static maps
    % Inputs:
    % - J    : optimization criterion [function handle]
    % - dt   : simulation step [s]
    % - N    : simulation lenght 
    % - f    : sinusoidal dither frequency [Hz]
    % - A    : sinusoidal dither amplitude [rad]
    % - fc   : HPF cut-off frequency [Hz]
    % - K    : integrator gain
    % - u0   : input initial value
    % - Ahess: hessian estimate dither amplitude
    % - fhess: hessian estimate dither frequency [HZ]
    % - wric : Ricatti filter parameter
    % - ddy0 : hessian inverse initial value
    % - AWA  : time variant AWA
    % Outputs:
    % - u      : control variable
    % - y      : criterion output
    % - dy     : criterion gradient estimate
    % - ddy    : hessian estimate
    % - ddy_inv: hessian inverse estimate
    
    global ship;

    % HPF
    bworder = 2;
    [b,a]   = butter(bworder, fc*dt, 'high');

    % Data structures
    u_hat   = [u0; zeros(N, 1)];
    u       = [u0; zeros(N, 1)];
    y       = zeros(N, 1);
    dy      = zeros(N, 1);
    hpf     = zeros(N, 1);
    ddy     = zeros(N, 1);
    ddy_inv = [ddy0; zeros(N, 1)]; % output of ricatti filter

    tic
    for i = 1:N
        if rem(i, 10) == 0
            fprintf("Iteration %d\n", i);
        end
        t = i*dt;

        ship.yaw = AWA(i);
        y(i) = J(u(i));
        
        % HPF
        if i >= bworder+1
            for j = 1:bworder+1
                hpf(i) = hpf(i) + b(j)*y(i-j+1);
            end
    
            for j = 2:bworder+1
                hpf(i) = hpf(i) - a(j)*hpf(i-j+1);
            end
            
            hpf(i) = 1/a(1) * hpf(i);
        end
        
        dy(i)        = hpf(i) * sin(2*pi*f*t); % M(t)y_hp
        
        ddy(i)       = hpf(i) * Ahess * (sin(2*pi*fhess*t)^2 - 0.5); % \hat{H} = N(t)y_hp
        
        ddy_inv(i+1) = ddy_inv(i) + dt * (wric * ddy_inv(i) - ...
            wric * ddy_inv(i)^2 * ddy(i)); % Euler discretization of Ricatti equation
    
        u_hat(i+1) = u_hat(i) - ddy_inv(i+1) * K * dt * dy(i); % single integrator
        u(i+1)     = u_hat(i+1) + A * sin(2*pi*f*t);
        
        if u(i+1) > pi || u(i+1) < -pi 
            break
        end
    end
    toc
end