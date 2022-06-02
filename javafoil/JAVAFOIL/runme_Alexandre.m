%--------------------------------------------------------------------------
% Performce optimization of thrust coefficeint CT(CL,CD,AWA) for any 
% number on identical profiles spaced in an arbitrary way on the hull centerline.
%
% The calculations are based on JavaFoil and includes semi-empirical
% estimation of turbulence transition and separation.
%
% Output is:
% CT - Thrust coefficent
% CL - for complete Rig setup.
% CD - For complete Rig setup
%-------------------------------------------------------------------------- 
clear; clc; close all;  addpath JavaFoil; addpath Foils;  
global ship counter;
fprintf('-------------------------------------------------------------\n');
figure(1);clf;figure(2);clf;

ship   = Ship(200); % (yaw,Pivot )Create the ship object
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord
R1     = Rig(26,0); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord,
% R1.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

R2     = Rig(75,0); % pivot x,y, 
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% R2.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord,

% R3     = Rig(122,0); % pivot x,y,  
% R3.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% % R3.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord
% 
% R4     = Rig( 170,0); % pivot x,y,  
% R4.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% % R4.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

ship.addRig(R1);  
ship.addRig(R2);  
% ship.addRig(R3);  
% ship.addRig(R4);

calc_scale();
ship.yaw = deg2rad(0);
% X = -ship.yaw*[1 1 1 1] + deg2rad(10) + [0 1 0 0]*deg2rad(20); % The trim angles
% [obj cl CT] = calc_objective(X);

%% Test section
ship.yaw = deg2rad(10);
[~,~,~] = calc_objective([deg2rad(0), deg2rad(0)]);

%% cL, cD curve - 1D
ship.yaw = deg2rad(0);
AoA = deg2rad(linspace(-20, 45, 66)); % sheeting angle = AoA
cL = zeros(length(AoA), 1); 
cD = zeros(length(AoA), 1); 
for i = 1:length(AoA)
    [~,cL(i),~,cD(i)] = calc_objective([AoA(i), 0, 0, 0]);
end

figure(10); clf(10);
title('cL(AoA)')
plot(rad2deg(AoA), cL, 'Linewidth', 2)
xlabel('AoA (º)');
ylabel('cL');

figure(11); clf(11); hold on;
title('cL(AoA), cD(AoA)')
plot(rad2deg(AoA), cL, 'Linewidth', 2)
plot(rad2deg(AoA), cD, 'r--', 'Linewidth', 1)
xlabel('AoA (º)');
legend('cL', 'cD');

%% cT curves - 1D
sheeting_angle = linspace(deg2rad(-180), deg2rad(180));
yaw = deg2rad([150, 120, 90, 60]);
cT = zeros(length(yaw), length(sheeting_angle)); % Thrust
cD = zeros(length(yaw), length(sheeting_angle)); % Drag

for i = 1:length(yaw)
    ship.yaw = yaw(i);
    tic
    for j = 1:length(sheeting_angle)
        [~,~, cT(i, j), cD(i, j)] = calc_objective(sheeting_angle(j));
    end
    toc
    figure(10+i); clf(10+i); hold on;
    title(['AWA = ', num2str(rad2deg(yaw(i))), 'º']);
    plot(rad2deg(sheeting_angle), cT(i, :), 'Linewidth', 1.8);
    plot(rad2deg(sheeting_angle), cD(i, :),'r--', 'Linewidth', 0.6);
    xlabel('Sheeting Angle (º)');
    ylabel('cT');
    legend('cT', 'cD')
end

figure(10+i+1); clf(10+i+1); hold on;
for i = 1: length(yaw)
    plot(rad2deg(sheeting_angle), cT(i, :), 'Linewidth', 2);
end
xlabel('Sheeting Angle (º)');
ylabel('cT');
legendCell = cellstr(num2str(rad2deg(yaw'), 'AWA=%.0dº'));
legend(legendCell);

%% cT surface 1D - Time-variant AWA
sheeting_angle = linspace(deg2rad(-180), deg2rad(180), 180);
yaw = linspace(deg2rad(-180), deg2rad(180), 180);
cT = zeros(length(yaw), length(sheeting_angle));

tic
for i = 1: length(yaw)
    ship.yaw = yaw(i);
    for j = 1:length(sheeting_angle)
        tic
        [~,~,cT(i, j),~] = calc_objective(sheeting_angle(j));
        t = toc;
    end
end
toc

figure(10); clf(10); hold on;
[X,Y] = meshgrid(rad2deg(yaw), rad2deg(sheeting_angle));

% Uncomment lines below to save struct
% data.cT = cT;
% data.x_grid = X';
% data.y_grid = Y';
% save('cT_SA_AWA.mat', 'data')

surf(X', Y', cT)
c = colorbar;
c.Label.String = 'cT';

title('cT(AWA, $\delta_s$)', 'Interpreter', 'Latex')
xlabel('AWA [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle $\delta_s$ [deg]', 'Interpreter', 'Latex');

%% cT surface 1D - Print surface
load('cT_SA_AWA.mat')

figure(10); clf(10); hold on;
surf(data.x_grid, data.y_grid, data.cT)
c = colorbar;
c.Label.String = 'cT';

for i = 1:size(data.y_grid, 1)
%     [~, locs] = findpeaks(data.cT(i, :), 'MinPeakDistance', 3.6*2)
    [~, locs] = max(data.cT(i, :));
    plot3(data.x_grid(i,locs), data.y_grid(i,locs), data.cT(i, locs), 'r.', 'Markersize', 15)
end

title('cT(AWA, $\delta_s$)', 'Interpreter', 'Latex')
xlabel('AWA [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle $\delta_s$ [deg]', 'Interpreter', 'Latex');

%% cT surface 2D - Fixed AWA
ship.yaw = deg2rad(45);
ship.addRig(R2)
scale = calc_scale();

sheeting_angle_1 = linspace(deg2rad(-180), deg2rad(180));
sheeting_angle_2 = linspace(deg2rad(-180), deg2rad(180));
cT = zeros(length(sheeting_angle_1), length(sheeting_angle_2));

tic
for i = 1: length(sheeting_angle_1)
    for j = 1:length(sheeting_angle_2)
        [~,~,cT(i, j),~] = calc_objective([sheeting_angle_1(i), sheeting_angle_2(j)]);
    end
end
toc

figure(10); clf(10); hold on;
[X,Y] = meshgrid(rad2deg(sheeting_angle_1), rad2deg(sheeting_angle_2));

% Uncomment lines below to save struct
% data.cT = cT;
% data.x_grid = X';
% data.y_grid = Y';
% save('cT1_cT2_SA.mat', 'data')

surf(X', Y', cT)
c = colorbar;
c.Label.String = 'cT';

% Hardcode AWA for title
title('cT($\delta_s^1, \delta_s^2$) | AWA = 45', 'Interpreter', 'Latex')
xlabel('sheeting angle 1 $\delta_s^1$ [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle 2 $\delta_s^2$ [deg]', 'Interpreter', 'Latex');

%% cT surface 2D - Fixed AWA - Print surface
load('cT1_cT2_SA.mat')

cT = data.cT;
cT(abs(cT) > 2.5) = nan;

figure(10); clf(10); hold on;
surf(data.x_grid, data.y_grid, cT)
c = colorbar;
c.Label.String = 'cT';

title('cT($\delta_s^1, \delta_s^2$) | AWA = 45', 'Interpreter', 'Latex')
xlabel('sheeting angle 1 $\delta_s^1$ [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle 2 $\delta_s^2$ [deg]', 'Interpreter', 'Latex');

%% cT data - 2 sail, time-variant AWA
sheeting_angle_1 = linspace(deg2rad(-180), deg2rad(180), 3);
sheeting_angle_2 = linspace(deg2rad(-180), deg2rad(180), 3);
AWA = linspace(deg2rad(-180), deg2rad(180), 2);
cT = zeros(length(AWA), length(sheeting_angle_1), length(sheeting_angle_2));

tic
for k = 1:length(AWA)
    ship.yaw = AWA(k);
    for i = 1: length(sheeting_angle_1)
        for j = 1:length(sheeting_angle_2)
            [~,~,cT(k, i, j),~] = calc_objective([sheeting_angle_1(i), sheeting_angle_2(j)]);
        end
    end
end
toc

% Uncomment lines below to save struct
data.cT = cT;
data.awa = AWA;
data.sheeting_angle_1 = sheeting_angle_1;
data.sheeting_angle_2 = sheeting_angle_2;
save('cT_2D_SA_AWA.mat', 'data')

%% Get slope reference
load('cT_SA_AWA.mat')
AWA = 45;
[~, i] = min(abs(data.x_grid(:, 1) - AWA));
cT_data = data.cT(i,  data.y_grid(1,:) < 20); % For AWA > 0 
[~, i] = max(cT_data);
data.y_grid(1, i);
slope = (cT_data(i-1) - cT_data(i)) / (data.y_grid(1, i-1) - data.y_grid(1, i));

%% References for 7M data - 1D
fig_cnt = 1;
dir = 'plots\7m_data_';

% TACKING AWA
sheeting_angle = linspace(deg2rad(-70), deg2rad(70), 70);
yaw = linspace(deg2rad(-80), deg2rad(80), 80);
cT = zeros(length(yaw), length(sheeting_angle));

tic
for i = 1: length(yaw)
    ship.yaw = yaw(i);
    for j = 1:length(sheeting_angle)
        [~,~,cT(i, j),~] = calc_objective(sheeting_angle(j));
    end
end
toc

figure(fig_cnt); clf(fig_cnt); hold on;
[X,Y] = meshgrid(rad2deg(yaw), rad2deg(sheeting_angle));

% Uncomment lines below to save struct
data.cT = cT;
data.x_grid = X';
data.y_grid = Y';
save(strcat(dir,'tacking\cT_SA_AWA.mat'), 'data');

surf(X', Y', cT)
c = colorbar;
c.Label.String = 'cT';

for i = 1:size(data.y_grid, 1)
    [~, locs] = max(data.cT(i, :));
    plot3(data.x_grid(i,locs), data.y_grid(i,locs), data.cT(i, locs), 'r.', 'Markersize', 15)
end

title('cT(AWA, $\delta_s$)', 'Interpreter', 'Latex')
xlabel('AWA [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle $\delta_s$ [deg]', 'Interpreter', 'Latex');
savefig(strcat(dir,'tacking\cT_SA_AWA.fig'))
hold off;
fig_cnt = fig_cnt + 1;

% 100 AWA
sheeting_angle = linspace(deg2rad(-110), deg2rad(-60), 50);
yaw = linspace(deg2rad(80), deg2rad(125), 45); 
cT = zeros(length(yaw), length(sheeting_angle));

tic
for i = 1: length(yaw)
    ship.yaw = yaw(i);
    for j = 1:length(sheeting_angle)
        [~,~,cT(i, j),~] = calc_objective(sheeting_angle(j));
    end
end
toc

figure(fig_cnt); clf(fig_cnt); hold on;
[X,Y] = meshgrid(rad2deg(yaw), rad2deg(sheeting_angle));

% Uncomment lines below to save struct
data.cT = cT;
data.x_grid = X';
data.y_grid = Y';
save(strcat(dir,'AWA_100\cT_SA_AWA.mat'), 'data');

surf(X', Y', cT)
c = colorbar;
c.Label.String = 'cT';

for i = 1:size(data.y_grid, 1)
    [~, locs] = max(data.cT(i, :));
    plot3(data.x_grid(i,locs), data.y_grid(i,locs), data.cT(i, locs), 'r.', 'Markersize', 15)
end

title('cT(AWA, $\delta_s$)', 'Interpreter', 'Latex')
xlabel('AWA [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle $\delta_s$ [deg]', 'Interpreter', 'Latex');
savefig(strcat(dir,'AWA_100\cT_SA_AWA.fig'))
hold off;
fig_cnt = fig_cnt + 1;

%% Smoothen reference surfance for 7M data - 1D
% Requires Image Processing Toolbox - To be installed
% AWA = +/- 45
load('plots\7m_data_tacking\cT_SA_AWA.mat')

surf(data.x_grid, data.y_grid, medfilt2(data.cT))
c = colorbar;
c.Label.String = 'cT';

% data.cT = medfilt2(data.cT);
% save('plots\7m_data_tacking\cT_SA_AWA.mat', 'data')

% AWA 100
load('plots\7m_data_AWA_100\cT_SA_AWA.mat')

surf(data.x_grid, data.y_grid, medfilt2(data.cT))
c = colorbar;
c.Label.String = 'cT';

% data.cT = medfilt2(data.cT);
% save('plots\7m_data_AWA_100\cT_SA_AWA.mat', 'data')

%% References for 7M data - 2D
dir = 'plots\7m_data_';

% TACKING AWA
sheeting_angle_1 = linspace(deg2rad(-90), deg2rad(90), 60); % res=3º
sheeting_angle_2 = linspace(deg2rad(-90), deg2rad(90), 60); % res=3º
yaw = linspace(deg2rad(-80), deg2rad(80), 40); % res=4º

cT = zeros(length(AWA), length(sheeting_angle_1), length(sheeting_angle_2));

tic
for k = 1:length(AWA)
    ship.yaw = AWA(k);
    for i = 1: length(sheeting_angle_1)
        for j = 1:length(sheeting_angle_2)
            [~,~,cT(k, i, j),~] = calc_objective([sheeting_angle_1(i), sheeting_angle_2(j)]);
        end
    end
end
toc

% Uncomment lines below to save struct
data.cT = cT;
data.AWA = AWA;
data.sheeting_angle_1 = sheeting_angle_1;
data.sheeting_angle_2 = sheeting_angle_2;
save(strcat(dir,'tacking\cT_2D_SA_AWA.mat'), 'data');

% 100 AWA.
sheeting_angle_1 = linspace(deg2rad(-125), deg2rad(-20), 50); % res=2º
sheeting_angle_2 = linspace(deg2rad(-125), deg2rad(-20), 50); % res=2º
AWA = linspace(deg2rad(80), deg2rad(125), 20); % res=2º

cT = zeros(length(AWA), length(sheeting_angle_1), length(sheeting_angle_2));

tic
for k = 1: length(AWA)
    ship.yaw = AWA(k);
    for i = 1: length(sheeting_angle_1)
        for j = 1:length(sheeting_angle_2)
            [~,~,cT(k, i, j),~] = calc_objective([sheeting_angle_1(i), sheeting_angle_2(j)]);
        end
    end
end
toc

% Uncomment lines below to save struct
data.cT = cT;
data.AWA = AWA;
data.sheeting_angle_1 = sheeting_angle_1;
data.sheeting_angle_2 = sheeting_angle_2;
save(strcat(dir,'AWA_100\cT_2D_SA_AWA.mat'), 'data');

%% Smoothen reference surfance for 7M data - 2D
% Requires Image Processing Toolbox - To be installed

% AWA = +/- 45
load('plots\7m_data_AWA_100\cT_2D_SA_AWA.mat')
data.cT = medfilt2(data.cT);
save('plots\7m_data_AWA_100\cT_2D_SA_AWA.mat', 'data')

% AWA 100
load('plots\7m_data_tacking\cT_2D_SA_AWA.mat')
data.cT = medfilt2(data.cT);
save('plots\7m_data_tacking\cT_2D_SA_AWA.mat', 'data')