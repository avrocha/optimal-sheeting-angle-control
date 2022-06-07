% Get 2D references w/ search-space constrained to real-data
% characteristics

clear; clc; close all;  
addpath JavaFoil; addpath Foils; addpath lib
global ship;
fprintf('-------------------------------------------------------------\n');
figure(1); clf; figure(2); clf;

ship   = Ship(200); % (yaw,Pivot) Create the ship object
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord
R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% R1.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

R2     = Rig(75,0); % pivot x,y, 
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% R2.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

ship.addRig(R1);  
ship.addRig(R2);

calc_scale();

tacking_data_dir  = 'data\measured_data\awa_pm_45\cT_2D.mat';
constant_data_dir = 'data\measured_data\awa_100\cT_2D.mat';

% 4-Core multiprocessing
parpool('local', 4);

% TACKING AWA
sheeting_angle_1 = linspace(deg2rad(-90), deg2rad(90), 60); % res=3º
sheeting_angle_2 = linspace(deg2rad(-90), deg2rad(90), 60); % res=3º
AWA = linspace(deg2rad(-80), deg2rad(80), 40); % res=4º

% Uncomment lines below to save data
data.AWA = AWA;
data.sheeting_angle_1 = sheeting_angle_1;
data.sheeting_angle_2 = sheeting_angle_2;

cT = zeros(length(AWA), length(sheeting_angle_1), length(sheeting_angle_2));

diary 'data\measured_data\awa_pm_45\cT_2D_diary.txt'

localShip = ship;
tic
for k = 1:1%length(AWA)
    localShip.yaw = AWA(k);
    fprintf("Iteration k = %d | AWA = %dº\n", k, rad2deg(localShip.yaw))
    parfor i = 1:4%length(sheeting_angle_1)
        aux = zeros(1, length(sheeting_angle_2));
        for j = 1:4%length(sheeting_angle_2)
            aux(j) = getfield(calc_objective_mod([sheeting_angle_1(i), sheeting_angle_2(j)], localShip), 'cT');
        end
        cT(k, i, :) = aux;
    end
    
    % Uncomment lines below to save data
    data.cT = cT;
    data.last_idxs.k = k;
    save(tacking_data_dir, 'data');
end
toc

diary off

% 100 AWA
sheeting_angle_1 = linspace(deg2rad(-125), deg2rad(-20), 50); % res=2º
sheeting_angle_2 = linspace(deg2rad(-125), deg2rad(-20), 50); % res=2º
AWA = linspace(deg2rad(80), deg2rad(125), 20); % res=2º

cT = zeros(length(AWA), length(sheeting_angle_1), length(sheeting_angle_2));

% Uncomment lines below to save data
data.AWA = AWA;
data.sheeting_angle_1 = sheeting_angle_1;
data.sheeting_angle_2 = sheeting_angle_2;

diary 'data\measured_data\awa_100\cT_2D_diary.txt'

localShip = ship;
tic
for k = 1:length(AWA)
    fprintf("Iteration k = %d | AWA = %dº\n", k, rad2deg(localShip.yaw))
    localShip.yaw = AWA(k);
    parfor i = 1:length(sheeting_angle_1)
        aux = zeros(1, length(sheeting_angle_2));
        for j = 1:length(sheeting_angle_2)
            aux(j) = getfield(calc_objective_mod([sheeting_angle_1(i), sheeting_angle_2(j)], localShip), 'cT');
        end
        cT(k, i, :) = aux;
    end
    
    % Uncomment lines below to save data
    data.cT = cT;
    data.last_idxs.k = k;
    save(constant_data_dir, 'data');
end
toc

poolobj = gcp('nocreate');
delete(poolobj);

diary off
