% Get 2D references w/ search-space constrained to real-data
% characteristics

clearvars -except F*; 
clc; close all;  
addpath JavaFoil; addpath Foils; addpath lib
global ship;
fprintf('-------------------------------------------------------------\n');
figure(1); clf; figure(2); clf;

ship   = Ship(200); % (yaw,Pivot) Create the ship object
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord
R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

R2     = Rig(75,0); % pivot x,y, 
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

ship.addRig(R1);  
ship.addRig(R2);

calc_scale();

tacking_data_dir  = 'data\measured_data\awa_pm_45\cT_2D.mat';
constant_data_dir = 'data\measured_data\awa_100\cT_2D.mat';

% 4-Core multiprocessing
parpool('local', 4);

% Function handle
fun = @(delta_1, delta_2, localShip, task_id) getfield(calc_objective_mod([delta_1, delta_2], localShip, task_id), 'cT');

% TACKING AWA
% ------------------------------------------------------------------------
sheeting_angle_1 = linspace(deg2rad(-90), deg2rad(90), 60); % res=3º (60)
sheeting_angle_2 = linspace(deg2rad(-90), deg2rad(90), 60); % res=3º (60)
AWA = linspace(deg2rad(-80), deg2rad(80), 40); % res=4º (40)

% Uncomment lines below to save data
data.AWA = AWA;
data.sheeting_angle_1 = sheeting_angle_1;
data.sheeting_angle_2 = sheeting_angle_2;

L1 = length(sheeting_angle_1);
L2 = length(sheeting_angle_2);
cT = zeros(length(AWA), L1, L2);

diary 'data\measured_data\awa_pm_45\cT_2D_diary.txt'

tic
localShip = ship;

for k = 1:length(AWA)
    localShip.yaw = AWA(k);
    fprintf("Iteration k = %d | AWA = %dº\n", k, rad2deg(localShip.yaw))

    parfor i = 1:L1
        delta_1 = sheeting_angle_1(i);
        for j = 1:L2
            cT(k, i, j) = fun(delta_1, sheeting_angle_2(j), localShip, getCurrentTask().ID);
        end
    end
    
    % Uncomment lines below to save data
    data.cT = cT;
    data.last_idxs.k = k;
    save(tacking_data_dir, 'data');
end
toc

diary off

% 100 AWA
% ------------------------------------------------------------------------
sheeting_angle_1 = linspace(deg2rad(-125), deg2rad(-20), 50); % res=2º (50)
sheeting_angle_2 = linspace(deg2rad(-125), deg2rad(-20), 50); % res=2º (50)
AWA = linspace(deg2rad(80), deg2rad(125), 20); % res=2º (20)

% Uncomment lines below to save data
data.AWA = AWA;
data.sheeting_angle_1 = sheeting_angle_1;
data.sheeting_angle_2 = sheeting_angle_2;

L1 = length(sheeting_angle_1);
L2 = length(sheeting_angle_2);
cT = zeros(length(AWA), L1, L2);

diary 'data\measured_data\awa_100\cT_2D_diary.txt'

tic
localShip = ship;

for k = 1:length(AWA)
    localShip.yaw = AWA(k);
    fprintf("Iteration k = %d | AWA = %dº\n", k, rad2deg(localShip.yaw))
    
    parfor i = 1:L1
        delta_1 = sheeting_angle_1(i);
        for j = 1:L2
            cT(k, i, j) = fun(delta_1, sheeting_angle_2(j), localShip, getCurrentTask().ID);
        end
    end
    
    % Uncomment lines below to save data
    data.cT = cT;
    data.last_idxs.k = k;
    save(constant_data_dir, 'data');
end
toc

diary off

poolobj = gcp('nocreate');
delete(poolobj);

%% Interpolation