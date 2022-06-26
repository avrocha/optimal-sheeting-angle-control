% Get 4D references w/ search-space constrained to real-data
% characteristics
%% Init
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

data_dir  = 'data\measured_data\awa_pm_45\cT_4D.mat';
diary 'data\measured_data\awa_pm_45\cT_4D_diary.txt'

% 4-Core multiprocessing
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool('local', 4);
end

% Function handle
fun = @(delta_1, delta_2, delta_3, delta_4, localShip, task_id) getfield(calc_objective_mod([delta_1, delta_2, delta_3, delta_4], localShip, task_id), 'cT');

% 45 AWA
% ------------------------------------------------------------------------
if isfile(data_dir)
    load(data_dir);
    L_awa = size(data.AWA, 2);
    L_sa           = size(data.sheeting_angle, 2);
    AWA = data.AWA;
    sheeting_angle = data.sheeting_angle;
end

if ~isfile(data_dir) || data.last_idxs.i == 0 
    fprintf("[Warning] Axes and cT structure init.\n")
    pause(5)
    
    % Constrain search space to the optima neighborhood
    optima_data = load('data\optima_calc\4D_optima.mat');
    
    AWA = linspace(deg2rad(15), deg2rad(80), 30); % (approx) res=2º (30)
    L_awa = size(AWA, 2);
    
    sheeting_angle = zeros(4, 10, L_awa); % res=2º (10)
    L_sa           = size(sheeting_angle, 2);
    
    for i = 1:length(AWA)
        [~, idx] = min(abs(AWA(i) - optima_data.data.AWA));
        for j = 1:4
            sheeting_angle(j, :, i) = linspace(optima_data.data.X(idx, j) - deg2rad(10), optima_data.data.X(idx, j) + deg2rad(10), L_sa);
        end
    end
    
    cT = zeros(L_awa, L_sa, L_sa, L_sa, L_sa);
    
    % Init structure and save axis and empty cT
    data.AWA            = AWA;
    data.sheeting_angle = sheeting_angle;
    data.cT             = cT;
    data.last_idxs.i    = 0;
    
    save(data_dir, 'data');
end

%% Confirmation
fprintf("ITERATION IS STARTING AT INDEX: *%d*\n\n", data.last_idxs.i+1);
pause(5);

fprintf("[2º Confirmation] ITERATION IS STARTING AT INDEX: *%d*\n\n", data.last_idxs.i+1);
pause(5);

%% Iteration
tic
localShip = ship;

for i = data.last_idxs.i+1:L_awa
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
    save(data_dir, 'data');
    
    if any(data.last_idxs.i == [10, 20])
        break;
    end
end
toc

%% Finish
diary off

poolobj = gcp('nocreate');
delete(poolobj);