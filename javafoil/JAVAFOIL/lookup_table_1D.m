% Get 2D references w/ search-space constrained to real-data
% characteristics

clear; clc; close all;  
addpath JavaFoil; addpath Foils; addpath lib
global ship;
fig_cnt = 1;
fprintf('-------------------------------------------------------------\n');
figure(1); clf; figure(2); clf;

ship   = Ship(200); % (yaw,Pivot) Create the ship object
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord
R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord

ship.addRig(R1);  

calc_scale();

tacking_data_dir  = 'data\measured_data\awa_pm_45\cT_1D.mat';
constant_data_dir = 'data\measured_data\awa_100\cT_1D.mat';

% 4-Core multiprocessing
% parpool('local', 4);

%% TACKING AWA
sheeting_angle = linspace(deg2rad(-80), deg2rad(80), 160); % res = 1º
AWA = linspace(deg2rad(-80), deg2rad(80), 160); % res = 1º

% Uncomment lines below to save data
data.AWA = AWA;
data.sheeting_angle = sheeting_angle;

cT = zeros(length(AWA), length(sheeting_angle));

diary 'data\measured_data\awa_pm_45\cT_1D_diary.txt'

localShip = ship;
tic
for k = 1:length(AWA)
    localShip.yaw = AWA(k);
    fprintf("Iteration k = %d | AWA = %dº\n", k, rad2deg(localShip.yaw))
    for i = 1:length(sheeting_angle)
        cT(k, i) = getfield(calc_objective_mod(sheeting_angle(i), localShip), 'cT');
    end
    % Uncomment lines below to save datad
    data.cT = cT;
    data.last_idxs.k = k;
    save(tacking_data_dir, 'data');
end
toc

diary off

% Print surface
figure(fig_cnt); clf(fig_cnt); hold on;
[X,Y] = meshgrid(rad2deg(AWA), rad2deg(sheeting_angle));
x_grid = X';
y_grid = Y';
surf(x_grid, y_grid, cT)
c = colorbar;
c.Label.String = 'cT';

for i = 1:size(y_grid, 1)
    [~, locs] = max(data.cT(i, :));
    plot3(x_grid(i,locs), y_grid(i,locs), data.cT(i, locs), 'r.', 'Markersize', 15)
end

title('cT(AWA, $\delta_s$)', 'Interpreter', 'Latex')
xlabel('AWA [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle $\delta_s$ [deg]', 'Interpreter', 'Latex');
savefig('data\measured_data\awa_pm_45\cT_1D_new.fig')
hold off;
fig_cnt = fig_cnt + 1;

%% 100 AWA
sheeting_angle = linspace(deg2rad(-120), deg2rad(-50), 70); % res 1º
AWA = linspace(deg2rad(80), deg2rad(125), 45); % res 1º

% Uncomment lines below to save data
data.AWA = AWA;
data.sheeting_angle = sheeting_angle;

cT = zeros(length(AWA), length(sheeting_angle));

diary 'data\measured_data\awa_100\cT_1D_diary_new.txt'

localShip = ship;
tic
for k = 1:length(AWA)
    localShip.yaw = AWA(k);
    fprintf("Iteration k = %d | AWA = %dº\n", k, rad2deg(localShip.yaw))
    for i = 1:length(sheeting_angle)
        cT(k, i) = getfield(calc_objective_mod(sheeting_angle(i), localShip), 'cT');
    end
    
    % Uncomment lines below to save data
    data.cT = cT;
    data.last_idxs.k = k;
    save(constant_data_dir, 'data');
end
toc

% poolobj = gcp('nocreate');
% delete(poolobj);

diary off

% Print surface
figure(fig_cnt); clf(fig_cnt); hold on;
[X,Y] = meshgrid(rad2deg(AWA), rad2deg(sheeting_angle));
x_grid = X';
y_grid = Y';
surf(x_grid, y_grid, cT)
c = colorbar;
c.Label.String = 'cT';

for i = 1:size(y_grid, 1)
    [~, locs] = max(data.cT(i, :));
    plot3(x_grid(i,locs), y_grid(i,locs), data.cT(i, locs), 'r.', 'Markersize', 15)
end

title('cT(AWA, $\delta_s$)', 'Interpreter', 'Latex')
xlabel('AWA [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle $\delta_s$ [deg]', 'Interpreter', 'Latex');
savefig('data\measured_data\awa_100\cT_1D_new.fig')
hold off;
fig_cnt = fig_cnt + 1;