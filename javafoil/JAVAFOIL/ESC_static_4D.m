%---
% Extremum seeking controller for a 1D static map
% J(\theta) = cT(sheeting_angle)
%---
clear; clc; close all;  
addpath JavaFoil;  addpath Foils;
global ship counter;
fprintf('-------------------------------------------------------------\n');

% Init configs
ship   = Ship(200); % (yaw,Pivot )Create the ship object
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord

% Rigs
R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord,

R2     = Rig(75,0); % pivot x,y, 
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

R3     = Rig(122,0); % pivot x,y,  
R3.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

R4     = Rig( 170,0); % pivot x,y,  
R4.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

ship.addRig(R1);  
ship.addRig(R2);  
ship.addRig(R3);  
ship.addRig(R4);
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
sheet_angle_0 = ones(4,1) * deg2rad(-40);
 
% Data structures
sheet_angle_hat = [sheet_angle_0; zeros(N, 4)];
sheet_angle     = [sheet_angle_0; zeros(N, 4)];
y               = zeros(N, 1);
grad_est        = zeros(N, 4);
filter_out      = zeros(N, 4);

ship.yaw = deg2rad(45);
tic
for n = 1:4
    for i = 1:N
        if rem(i, 10) == 0
            fprintf("Iteration %d\n", i);
        end
        t = i*dt;
        y(i*n) = J(sheet_angle(i, n));

        if i >= bworder+1
            for j = 1:bworder+1
                filter_out(i, n) = filter_out(i, n) + b(j)*y((i-j+1)*n);
            end

            for j = 2:bworder+1
                filter_out(i, n) = filter_out(i, n) - a(j)*filter_out(i-j+1, n);
            end

            filter_out(i, n) = 1/a(1) * filter_out(i, n);
        end

        grad_est(i, n) = filter_out(i, n) * sin(2*pi*f*t);

        sheet_angle_hat(i+1, n) = sheet_angle_hat(i, n) + K*dt*grad_est(i, n); % single integrator

        sheet_angle(i+1, n) = sheet_angle_hat(i+1, n) + A*sin(2*pi*f*t);

        if sheet_angle(i+1, n) > pi || sheet_angle(i+1, n) < -pi 
            break
        end
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