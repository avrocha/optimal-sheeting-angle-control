%---
% Extremum seeking controller for a 1D static map
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

% ESC params
% Simulation
fs = 100; % sampling frequency (Hz)
dt = 1/fs;
T  = 5;
N = length(0:dt:T);

% Perturbation
f = 10; % perturbation frequency (Hz)
A = deg2rad(2); % perturbation amplitude (rad)

% HPF
fc      = 5; % cutoff frequency (Hz) (has to be high enough)
bworder = 2;
[b,a]   = butter(bworder, fc*dt, 'high');

% Integrator
K             = 5; % gain (>0 since extremum is maximum)
sheet_angle_0 = deg2rad(-15);
 
% Data structures
sheet_angle_hat = [sheet_angle_0; zeros(N, 1)];
sheet_angle     = [sheet_angle_0; zeros(N, 1)];
y               = zeros(N, 1);
grad_est        = zeros(N, 1);
filter_out      = zeros(N, 1);

ship.yaw = deg2rad(45);
tic
for i = 1:N
    if rem(i, 10) == 0
        fprintf("Iteration %d\n", i);
    end
    t = i*dt;
    y(i) = J(sheet_angle(i));

    if i >= bworder+1
        for j = 1:bworder+1
            filter_out(i) = filter_out(i) + b(j)*y(i-j+1);
        end

        for j = 2:bworder+1
            filter_out(i) = filter_out(i) - a(j)*filter_out(i-j+1);
        end
        
        filter_out(i) = 1/a(1) * filter_out(i);
    end
    
    grad_est(i) = filter_out(i) * sin(2*pi*f*t);
    
    sheet_angle_hat(i+1) = sheet_angle_hat(i) + K*dt*grad_est(i); % single integrator
        
    sheet_angle(i+1) = sheet_angle_hat(i+1) + A*sin(2*pi*f*t);
    
    if sheet_angle(i+1) > pi || sheet_angle(i+1) < -pi 
        break
    end
end
toc

%% Plots
figure(1); clf(1); hold on;
plot(0:dt:T, rad2deg(sheet_angle(1:N)), 'b-', 'Linewidth', 2)
plot(0:dt:T, -25*ones(N,1), 'r--', 'Linewidth', 1)
xlabel('t (s)'), ylabel('$\delta_s$', 'Interpreter', 'Latex')

figure(2); clf(2); hold on;
plot(0:dt:T, y, 'b-', 'Linewidth', 2)
xlabel('t (s)'), ylabel('$cT$', 'Interpreter', 'Latex')

figure(3); clf(3); hold on;
plot(0:dt:T, grad_est, 'b-', 'Linewidth', 2)
xlabel('t (s)'), ylabel('$\hat{\zeta}$', 'Interpreter', 'Latex')

figure(4); clf(4); hold on;
plot(0:dt:T, filter_out, 'b-', 'Linewidth', 2)
xlabel('t (s)'), ylabel('Filter Output')