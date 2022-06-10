%---
% [Frequency selection script implemented in parallel - REQUIRES PARALLEL LIBRARY]
%---
% Copyright: Alexandre Vieira da Rocha

clear; clc; close all;  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

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

% Simulation
fs = 2; % sampling frequency (Hz)
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
sheet_angle_0 = [deg2rad(-35)*ones(1, 2);
                 deg2rad(-80)*ones(1, 2);
                 deg2rad(-125)*ones(1, 2)];

% Accumulated cT (GB | NB)
cT_cum = zeros(size(AWA, 1), size(f, 1), 2);

p = gcp('nocreate');
if isempty(p)
    p = parpool('local', 4);
end

localShip = ship;

% GB 
J = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 2), 'cT')); % task_id = 2

gb_fun = @(f, A, fc_hp, K, x0, AWA) gbesc(localShip, J, dt, N, f, A, fc_hp, fc_hp, K, x0, AWA, false);
loop   = @() gb_loop(gb_fun, AWA, f, sheet_angle_0);
F1     = parfeval(loop, 2);

% NB 
J = @(sheeting_angle, ship)(getfield(calc_objective_mod(sheeting_angle, ship, 3), 'cT')); % task_id = 3

nb_fun = @(f, A, fc_hp, K, x0, wric, ric_0, AWA) nbesc(localShip, J, dt, N, f, A, fc_hp, fc_hp, K, x0, wric, ric_0, AWA, false);
loop    = @() nb_loop(nb_fun, AWA, f, sheet_angle_0);
F2      = parfeval(loop, 2);

F1.State
F2.State

%% Results
% load("data/cT_cum_freq_selection.mat");

[gb.cT, gb.cT_cum] = fetchOutputs(F1);
[nb.cT, nb.cT_cum] = fetchOutputs(F2);

poolobj = gcp('nocreate');
delete(poolobj);

diary 'data\freq_selection\diary.txt'

fprintf("GB Results - Frequency Selection\n")
for k = 1:size(AWA, 1)
    for i = 1:size(f, 1)
        fprintf("\t AWA = %s, f = %s, cT_cum = %f\n", num2str(AWA(k, 1)), num2str(f(i, :)), gb.cT_cum(k,i))
    end
end

fprintf("\n--------------||--------------\n\n")

fprintf("NB Results - Frequency Selection\n")
for k = 1:size(AWA, 1)
    for i = 1:size(f, 1)
        fprintf("\t AWA = %s, f = %s, cT_cum = %f\n", num2str(AWA(k, 1)), num2str(f(i, :)), nb.cT_cum(k,i))
    end
end

diary off

save("data\freq_selection\cT_cum_freq_selec.mat", "gb", "nb");

%% Loops (1 per worker)

% GB 
function [cT, cT_cum] = gb_loop(fun, AWA, f, sheet_angle_0)
    % Fixed params
    A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
    K             = diag([0.0750, 0.0750]); % gain (>0 since extremum is maximum)

    cT_cum = zeros(size(AWA, 1), size(f, 1));

    for k = 1:size(AWA, 1)
        for i = 1:size(f, 1)
            fc_hp = 0.7 * min(f(i, :)); % HPF cutoff freq

            [~, cT, ~] = fun(f(i, :)', A, fc_hp, K, sheet_angle_0(k, :)', AWA(k, :));
            
            cT_cum(k, i) = sum(cT, 2);
        end
    end
end

% NB
function [cT, cT_cum] = nb_loop(fun, AWA, f, sheet_angle_0)
    % Fixed params
    A             = [deg2rad(2); deg2rad(2)]; % dither amplitude
    K             = diag([0.0025, 0.0025]); % gain (>0 since extremum is maximum)
    wric          = 2 * pi * 0.003; % ricatti filter parameter: 0.05f (0.01f) without (with) LPF
    ric_0         = diag([-30, -30]);
    
    cT_cum = zeros(size(AWA, 1), size(f, 1));
    
    for k = 1:size(AWA, 1)
        for i = 1:size(f, 1)
            fc_hp = 0.7 * min(f(i, :)); % HPF cutoff freq
    
            [~, cT, ~, ~, ~] = fun(f(i, :)', A, fc_hp, K, sheet_angle_0(k, :)', wric, ric_0, AWA(k, :));

            cT_cum(k, i)  = sum(cT, 2);
        end
    end
end
