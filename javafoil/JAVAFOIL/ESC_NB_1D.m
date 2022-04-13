%---
% Newton-based extremum seeking controller for a 1D static map
% J(\theta) = cT(sheeting_angle)
%---
clear; clc; close all;  
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

% Simulation
fs = 1; % sampling frequency (Hz)
dt = 1/fs;
T  = 300;
N = length(0:dt:T);

% Perturbation
f = 0.1; % perturbation frequency (Hz)
A = deg2rad(2); % perturbation amplitude (rad)

% HPF
fc      = 0.05; % cutoff frequency (Hz) (has to be high enough)
bworder = 2;
[b,a]   = butter(bworder, fc*dt, 'high');

% Integrator
K             = 0.0025; % gain (>0 since extremum is maximum)
sheet_angle_0 = deg2rad(-15);
 
% Hessian estimate N(t)
Ahess = 16 / (rad2deg(A)^2);
fhess = f;

% Ricatti filter
wric = 2 * pi * 0.05 * f;

% Data structures
sheet_angle_hat = [sheet_angle_0; zeros(N, 1)];
sheet_angle     = [sheet_angle_0; zeros(N, 1)];
y               = zeros(N, 1);
grad_est        = zeros(N, 1);
filter_out      = zeros(N, 1);
hessian_est     = zeros(N, 1);
% [WIP] Test different inv(hess) init values
hessian_inv_est = [-100; zeros(N, 1)]; % output of ricatti filter

ship.yaw = deg2rad(45); % AWA
tic
for i = 1:N
    if rem(i, 10) == 0
        fprintf("Iteration %d\n", i);
    end
    t = i*dt;
    y(i) = J(sheet_angle(i));
    
    % HPF
    if i >= bworder+1
        for j = 1:bworder+1
            filter_out(i) = filter_out(i) + b(j)*y(i-j+1);
        end

        for j = 2:bworder+1
            filter_out(i) = filter_out(i) - a(j)*filter_out(i-j+1);
        end
        
        filter_out(i) = 1/a(1) * filter_out(i);
    end
    
    grad_est(i)          = filter_out(i) * sin(2*pi*f*t); % M(t)y_hp
    
    hessian_est(i)       = filter_out(i) * Ahess * (sin(2*pi*fhess*t)^2 - 0.5); % \hat{H} = N(t)y_hp
    
    hessian_inv_est(i+1) = hessian_inv_est(i) + dt * (wric * hessian_inv_est(i) - ...
        wric * hessian_inv_est(i)^2 * hessian_est(i)); % Euler discretization of Ricatti equation

    sheet_angle_hat(i+1) = sheet_angle_hat(i) - hessian_inv_est(i+1) * K * dt * grad_est(i); % single integrator
    sheet_angle(i+1)     = sheet_angle_hat(i+1) + A*sin(2*pi*f*t);
    
    if sheet_angle(i+1) > pi || sheet_angle(i+1) < -pi 
        break
    end
end
toc

%% Plots
% % Directory to save plots
% filename = 'plots\NB_ESC\';
% 
% % Save last params in diary.txt
% fileID = fopen(strcat(filename,'diary.txt'),'a');
% fprintf(fileID, "----%s----\n", datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
% fprintf(fileID,['Params:\n'...
%                 'HPF: fc = %f\n'...
%                 'Main dither: A = %f, f = %d\n'...
%                 'Ricatti Filter: wric = %f\n'...
%                 'Hessian estimate dither: Ahess =%f, whess = %f\n'...
%                 'Integrator gain: K = %f\n'...
%                 'Initial input value: %f\n\n'], [fc; A; f; wric; Ahess; fhess; K; sheet_angle_0]);
% fclose(fileID);

figure(1); clf(1); hold on;
title('NB-ESC | Sheeting Angle')
plot(0:dt:T, rad2deg(sheet_angle(1:N)), 'b-', 'Linewidth', 2)
plot(0:dt:T, -25*ones(N,1), 'r--', 'Linewidth', 1)
xlabel('t (s)'), ylabel('$\delta_s$', 'Interpreter', 'Latex')
print('plots\1c.eps','-depsc'); 

figure(2); clf(2); hold on;
title('NB-ESC | Thrust Coeff')
plot(0:dt:T, y, 'b-', 'Linewidth', 2)
xlabel('t (s)'), ylabel('$cT$', 'Interpreter', 'Latex')

figure(3); clf(3); hold on;
title('NB-ESC | Gradient Estimate')
plot(0:dt:T, grad_est, 'b-', 'Linewidth', 2)
plot(0:dt:T, movmean(grad_est, 10), 'r--', 'Linewidth', 1.3)
xlabel('t (s)'), ylabel('$\hat{\zeta}$', 'Interpreter', 'Latex')

figure(4); clf(4); hold on;
title('NB-ESC | Hessian Estimate')
plot(0:dt:T, hessian_est, 'b-', 'Linewidth', 2)
plot(0:dt:T, movmean(hessian_est, 10), 'r--', 'Linewidth', 1.3)
xlabel('t (s)'), ylabel('$\hat{H}$', 'Interpreter', 'Latex')

figure(5); clf(5); hold on;
title('NB-ESC | Hessian Inverse Estimate')
plot(0:dt:T+1, hessian_inv_est, 'b-', 'Linewidth', 2)
xlabel('t (s)'), ylabel('$\hat{H}^{-1}$', 'Interpreter', 'Latex')