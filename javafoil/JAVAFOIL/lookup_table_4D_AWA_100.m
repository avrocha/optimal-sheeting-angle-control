% Get 4D references w/ search-space constrained to real-data
% characteristics

clearvars -except F*; 
clc; close all;  
addpath JavaFoil; addpath Foils; addpath lib
global ship;
fprintf('-------------------------------------------------------------\n');
figure(1); clf; figure(2); clf;

ship   = Ship(200); % (yaw,Pivot )Create the ship object
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord
R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord,

R2     = Rig(75,0); % pivot x,y, 
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

R3     = Rig(122,0); % pivot x,y,  
R3.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

R4     = Rig(170,0); % pivot x,y,  
R4.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

ship.addRig(R1);
ship.addRig(R2);
ship.addRig(R3);
ship.addRig(R4);
calc_scale();

tacking_data_dir  = 'data\measured_data\awa_pm_45\cT_4D.mat';
constant_data_dir = 'data\measured_data\awa_100\cT_4D.mat';

% 4-Core multiprocessing
parpool('local', 4);

% Function handle
fun = @(delta_1, delta_2, delta_3, delta_4, localShip, task_id) getfield(calc_objective_mod([delta_1, delta_2, delta_3, delta_4], localShip, task_id), 'cT');

% 100 AWA
% ------------------------------------------------------------------------
% Constrain search space to the optima neighborhood
optima_data = load('data\optima_calc\4D_optima.mat');

AWA = linspace(deg2rad(85), deg2rad(125), 20); % res=2º (20)
L_awa = size(AWA, 2);

sheeting_angle = zeros(4, 10, 20); % res=2º (10)
L_sa           = size(sheeting_angle, 2);

for i = 1:length(AWA)
    [~, idx] = min(abs(AWA(i) - optima_data.data.AWA));
    for j = 1:4
        sheeting_angle(j, :, i) = linspace(optima_data.data.X(idx, j) - deg2rad(10), optima_data.data.X(idx, j) + deg2rad(10), 10);
    end
end

% Uncomment lines below to save data
data.AWA = AWA;
data.sheeting_angle = sheeting_angle;

cT    = zeros(L_awa, L_sa, L_sa, L_sa, L_sa);

diary 'data\measured_data\awa_100\cT_4D_diary.txt'

tic
localShip = ship;

for i = 1:L_awa
    localShip.yaw = AWA(i);
    fprintf("Iteration k = %d | AWA = %dº\n", i, rad2deg(localShip.yaw))
    
    delta = sheeting_angle(:, :, i);
    parfor j = 1:L_sa
        for k = 1:L_sa
            for l = 1:L_sa
                for m = 1:L_sa
                    cT(i, j, k, l, m) = fun(delta(1, j), delta(2, k), delta(3, l), delta(4, m), localShip, getCurrentTask().ID);
                end
            end
        end
    end
    
    % Uncomment lines below to save data
    data.cT = cT;
    data.last_idxs.i = i;
    save(constant_data_dir, 'data');
end
toc

diary off

poolobj = gcp('nocreate');
delete(poolobj);

%% TACKING AWA [WIP]
% ------------------------------------------------------------------------
% sheeting_angle_1 = linspace(deg2rad(-90), deg2rad(90), 60); % res=3º (60)
% sheeting_angle_2 = linspace(deg2rad(-90), deg2rad(90), 60); % res=3º (60)
% AWA = linspace(deg2rad(-80), deg2rad(80), 40); % res=4º (40)
% 
% % Uncomment lines below to save data
% data.AWA = AWA;
% data.sheeting_angle_1 = sheeting_angle_1;
% data.sheeting_angle_2 = sheeting_angle_2;
% 
% L1 = length(sheeting_angle_1);
% L2 = length(sheeting_angle_2);
% cT = zeros(length(AWA), L1, L2);
% 
% diary 'data\measured_data\awa_pm_45\cT_4D_diary.txt'
% 
% tic
% localShip = ship;
% 
% for k = 1:length(AWA)
%     localShip.yaw = AWA(k);
%     fprintf("Iteration k = %d | AWA = %dº\n", k, rad2deg(localShip.yaw))
% 
%     parfor i = 1:L1
%         delta_1 = sheeting_angle_1(i);
%         for j = 1:L2
%             cT(k, i, j) = fun(delta_1, sheeting_angle_2(j), localShip, getCurrentTask().ID);
%         end
%     end
%     
%     % Uncomment lines below to save data
%     data.cT = cT;
%     data.last_idxs.k = k;
%     save(tacking_data_dir, 'data');
% end
% toc
% 
% diary off