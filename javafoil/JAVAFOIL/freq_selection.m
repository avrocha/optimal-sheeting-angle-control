%---
% [SCRIPT FOR FREQUENCY ORDER SELECTION]
% Gradient- and Newton-based extremum seeking controller for multivariate time-variant AWA
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
R2     = Rig(75,0); % pivot x,y,  
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile, x, y, dx, chord

ship.addRig(R1);
ship.addRig(R2);

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
T  = 250;
N  = length(0:dt:T);
    
% Data structures
% Periods
T = [5, 10, 20];
f = [1/T(1), 1/T(2);
     1/T(2), 1/T(1);
     1/T(2), 1/T(3);
     1/T(3), 1/T(2);
     1/T(1), 1/T(3);
     1/T(3), 1/T(1)];

% AWA
AWA = [deg2rad(45)*ones(1, N);
       deg2rad(90)*ones(1, N);
       deg2rad(135)*ones(1, N)];

% u0
sheet_angle_0 = [deg2rad(-40)*ones(1, 2);
                 deg2rad(-85)*ones(1, 2);
                 deg2rad(-130)*ones(1, 2)];

% Accumulated cT (GB | NB)
cT_cum = zeros(size(AWA, 1), size(f, 1), 2);

%% GB-ESC loop
fprintf("Gradient-based ESC selected\n.")

% fixed params
A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
K             = diag([0.0750, 0.0750]); % gain (>0 since extremum is maximum)

for k = 1:size(AWA, 1)
    for i = 1:size(f, 1)
        fc_hp = 0.9 * min(f(i, :)); % HPF cutoff freq
        [~, cT, ~] = gbesc(J, dt, N, f(i, :)', A, fc_hp, fc_hp, K, sheet_angle_0(k, :)', AWA(k, :), false);
        cT_cum(k, i, 1) = sum(cT, 2);
    end
end

%% NB-ESC loop
fprintf("Newton-based ESC selected\n.")

% Fixed params
A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
K             = diag([0.0025, 0.0025]); % gain (>0 since extremum is maximum)
wric          = 2 * pi * 0.05 * min(f); % ricatti filter parameter: 0.05f (0.01f) without (with) LPF
ric_0         = diag([-30, -30]);

for k = 1:size(AWA, 1)
    for i = 1:size(f, 1)
        fc_hp = 0.9 * min(f(i, :)); % HPF cutoff freq
        [~, cT, ~, ~, ~] = nbesc(J, dt, N, f(i, :)', A, fc_hp, fc_hp, K, sheet_angle_0(k, :)', wric, ric_0, AWA(k, :), false);
        cT_cum(k, i, 2) = sum(cT, 2);
    end
end

%% Results
fprintf("GB Results - Frequency Selection\n")
for k = 1:size(AWA, 1)
    for i = 1:size(f, 1)
        fprintf("\t AWA = %s, f = %s, cT_cum = %f\n", num2str(AWA(k, 1)), num2str(f(i, :)), cT_cum(k,i,1))
    end
end

fprintf("\n--------------||--------------\n\n")

fprintf("NB Results - Frequency Selection\n")
for k = 1:size(AWA, 1)
    for i = 1:size(f, 1)
        fprintf("\t AWA = %s, f = %s, cT_cum = %f\n", num2str(AWA(k, 1)), num2str(f(i, :)), cT_cum(k,i,2))
    end
end


%% Functions
function [u, y, dy] = gbesc(J, dt, N, f, A, fc_hp, fc_lp, K, u0, AWA, lp_bool)
    % Inputs:
    % - J      : optimization criterion [function handle]
    % - dt     : simulation step [s]
    % - N      : simulation lenght 
    % - f      : sinusoidal dithers frequency [Hz]
    % - A      : sinusoidal dithers amplitude [rad]
    % - fc_hp  : HPF cut-off frequency [Hz]
    % - fc_lp  : LPPF cut-off frequency [Hz]
    % - K      : integrator gain
    % - u0     : initial sheeting angle [rad]
    % - AWA    : time variant AWA
    % - lp_bool: use LPF [boolean] 
    % Outputs:
    % - u : control variable
    % - y : criterion output
    % - dy: criterion gradient estimate
    
    global ship;

    bworder = 5;
    % HPF
    [bh,ah]   = butter(bworder, fc_hp*dt, 'high');
    % LPF
    [bl,al]   = butter(bworder, fc_lp*dt, 'low');

    % Data structures
    n     = length(f);
    u_hat = [u0, zeros(n, N)];
    u     = [u0, zeros(n, N)];
    y     = zeros(1, N);
    dy    = zeros(n, N);
    hpf   = zeros(1, N);
    zeta  = zeros(n, N);
    lpf   = zeros(n, N);
    
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
   
        ship.yaw = AWA(i);
        
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
            zeta(:, i) = hpf(i) * sin(2*pi*f*t);
            
            % LPF
            if i >= bworder+1
                for j = 1:bworder+1
                    lpf(:, i) = lpf(:, i) + bl(j) .* zeta(:, i-j+1);
                end
        
                for j = 2:bworder+1
                    lpf(:, i) = lpf(:, i) - al(j) .* lpf(:, i-j+1);
                end
                
                lpf(:, i) = 1/al(1) .* lpf(:, i);
            else
                for j = 1:i
                    lpf(:, i) = lpf(:, i) + bl(j) .* zeta(:, i-j+1);
                end
        
                for j = 2:i
                    lpf(:, i) = lpf(:, i) - al(j) .* lpf(:, i-j+1);
                end
                
                lpf(:, i) = 1/al(1) .* lpf(:, i);
            end
    
            dy(:, i) = lpf(:, i);
        
        else
            dy(:, i) = hpf(i) * sin(2*pi*f*t);

        end
        
        % Parameter estimate - Single Integrator
        u_hat(:, i+1) = u_hat(:, i) + dt * K * dy(:, i); % single integrator
        
        % Add dither
        u(:, i+1)     = u_hat(:, i+1) + A .* sin(2*pi*f*t);
        
        % Error condition
        if any(u(i+1) > pi) || any(u(i+1) < -pi) 
            break
        end
    end
    toc
end

function [u, y, dy, ddy, ddy_inv] = nbesc(J, dt, N, f, A, fc_hp, fc_lp, K, u0, wric, ddy0, AWA, lp_bool)
    % Inputs:
    % - J    : optimization criterion [function handle]
    % - dt   : simulation step [s]
    % - N    : simulation lenght 
    % - f    : sinusoidal dither frequency [Hz]
    % - A    : sinusoidal dither amplitude [rad]
    % - fc_hp  : HPF cut-off frequency [Hz]
    % - fc_lp  : LPF cut-off frequency [Hz]
    % - K    : integrator gain
    % - u0     : initial sheeting angle [rad]
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

    % Data structures
    n              = length(f);
    u_hat          = [u0, zeros(n, N)];
    u              = [u0, zeros(n, N)];
    y              = zeros(1, N);
    dy             = zeros(n, N);
    hpf            = zeros(1, N);
    zeta           = zeros(n, N);
    lpf_grad       = zeros(n, N);
    sigma          = zeros(n, n, N);
    lpf_hessian    = zeros(n, n, N);
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
   
        ship.yaw = AWA(i);
        
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
            % Gradient
            zeta(:, i) = hpf(i) .* sin(2*pi*f*t);
            
            % LPF - Gradient
            if i >= bworder+1
                for j = 1:bworder+1
                    lpf_grad(:, i) = lpf_grad(:, i) + bl(j) .* zeta(:, i-j+1);
                end
        
                for j = 2:bworder+1
                    lpf_grad(:, i) = lpf_grad(:, i) - al(j) .* lpf_grad(:, i-j+1);
                end
                
                lpf_grad(:, i) = 1/al(1) .* lpf_grad(:, i);
            else
                for j = 1:i
                    lpf_grad(:, i) = lpf_grad(:, i) + bl(j) .* zeta(:, i-j+1);
                end
        
                for j = 2:i
                    lpf_grad(:, i) = lpf_grad(:, i) - al(j) .* lpf_grad(:, i-j+1);
                end
                
                lpf_grad(:, i) = 1/al(1) .* lpf_grad(:, i);
            end
            
            dy(:, i) = lpf_grad(:, i);
            
            % Hessian 
            for j = 1:n
                for k = 1:n
                    sigma(j,k,i) = hpf(i) * Nhess{j,k}(t); % \hat{H} = N(t)y_hp
                end
            end

            % LPF - Hessian
            if i >= bworder+1
                for j = 1:bworder+1
                    lpf_hessian(:, :, i) = lpf_hessian(:, :, i) + bl(j).*sigma(:, :, i-j+1);
                end
        
                for j = 2:bworder+1
                    lpf_hessian(:, :, i) = lpf_hessian(:, :, i) - al(j).*lpf_hessian(:, :, i-j+1);
                end
                
                lpf_hessian(:, :, i) = 1/al(1) * lpf_hessian(:, :, i);
            else
                for j = 1:i
                    lpf_hessian(:, :, i) = lpf_hessian(:, :, i) + bl(j).*sigma(:, :, i-j+1);
                end
        
                for j = 2:i
                    lpf_hessian(:, :, i) = lpf_hessian(:, :, i) - al(j)*lpf_hessian(:, :, i-j+1);
                end
                
                lpf_hessian(:, :, i) = 1/al(1) .* lpf_hessian(:, :, i);
            end
    
            ddy(:, :, i) = lpf_hessian(:, :, i);
        
        else
            % Gradient
            dy(:, i) = hpf(i) * sin(2*pi*f*t);
            
            % Hessian
            for j = 1:n
                for k = 1:n
                    ddy(j,k,i) = hpf(i) * Nhess{j,k}(t); % \hat{H} = N(t)y_hp
                end
            end            
        end
        
        % Euler discretization of Ricatti equation    
        ddy_inv(:, :, i+1) = ddy_inv(:, :, i) + dt * (wric * ddy_inv(:, :, i) - ...
            wric * ddy_inv(:, :, i) * ddy(:, :, i) * ddy_inv(:, :, i)); 

        % Parameter estimate - Single Integrator
        u_hat(:, i+1) = u_hat(:, i) - dt * K * ddy_inv(:, :, i+1) * dy(:, i); 
        
        % Add dither
        u(:, i+1) = u_hat(:, i+1) + A .* sin(2*pi*f*t);
        
        % Error condition
        if any(u(:, i+1) > pi) || any(u(:, i+1) < -pi) 
            break
        end
    end
    toc
end