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
T  = 2;
N = length(0:dt:T);

% Perturbation
f = 25; % perturbation frequency (Hz)
A = deg2rad(1); % perturbation amplitude (rad)

% HPF
fc  = 0.02; % cutoff frequency (Hz)
bworder = 2;
[b,a] = butter(bworder, fc*dt, 'high');

% Integrator
K             = 10; % gain (>0 since extremum is maximum)
sheet_angle_0 = deg2rad(-30);

% Data structures
yaws = [45]; %, 90, 135];
sheet_angle = [sheet_angle_0; zeros(N, 1)];
y           = zeros(N, 1);
grad_est    = zeros(N, 1);

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
            grad_est(i) = grad_est(i) + b(j)*y(i-j+1);
        end

        for j = 2:bworder+1
            grad_est(i) = grad_est(i) - a(j)*grad_est(i-j+1);
        end
        
        grad_est(i) = 1/a(1) * grad_est(i);
    end
    
    grad_est(i) = grad_est(i) * sin(2*pi*f*t);
    
    sheet_angle(i+1) = sheet_angle(i) + K*dt*grad_est(i); % single integrator
    
    sheet_angle(i+1) = sheet_angle(i+1) + A*sin(2*pi*f*t);
end

%% Plots
figure(1); clf(1); hold on;
plot(0:dt:T, rad2deg(sheet_angle(1:N)), 'b-', 'Linewidth', 2)
plot(0:dt:T, -25*ones(N,1), 'r--', 'Linewidth', 1)
xlabel('t (s)'), ylabel('$\delta_s$', 'Interpreter', 'Latex')
ylim([-90,90])
figure(2); clf(2); hold on;
plot(0:dt:T, y, 'b-', 'Linewidth', 2)
xlabel('t (s)'), ylabel('$cT$', 'Interpreter', 'Latex')

figure(3); clf(3); hold on;
plot(0:dt:T, grad_est, 'b-', 'Linewidth', 2)
xlabel('t (s)'), ylabel('$\hat{\epsilon}$', 'Interpreter', 'Latex')