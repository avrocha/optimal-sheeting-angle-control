%--------------------------------------------------------------------------
% Performce optimization of thrust coefficeint CT(CL,CD,AWA) for any 
% number on identical profiles spaced in an arbitrary way on the hull centerline.
%
% The calculations are based on JavaFoil and includes semi-empirical
% estimation of turbulence transition and separation.
%
% Output is:
% CT - Thrust coefficent
% CL - for complete Rig setup
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

R3     = Rig(122,0); % pivot x,y,  
R3.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% R3.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

R4     = Rig( 170,0); % pivot x,y,  
R4.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% R4.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

ship.addRig(R1);  
% ship.addRig(R2);  
% ship.addRig(R3);  
% ship.addRig(R4);

scale = calc_scale();
ship.yaw = deg2rad(0);
% X = -ship.yaw*[1 1 1 1] + deg2rad(10) + [0 1 0 0]*deg2rad(20); % The trim angles
% [obj cl CT] = calc_objective(X);

% You can now choose to focus on ship cL (upwards in the picture) or cT
% which is the force in the ship direction.

%% Test section
ship.yaw = deg2rad(60);
[~,~,~] = calc_objective([deg2rad(-180), 0, 0, 0]);

%% cL, cD curve - 2D
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

%% cT curves - 2D
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

%% cT surface - 1 sail [WIP]
sheeting_angle = linspace(deg2rad(-180), deg2rad(180), 3);
yaw = linspace(deg2rad(0), deg2rad(180), 3);
cT = zeros(length(yaw), length(sheeting_angle));

tic
for i = 1: length(yaw)
    ship.yaw = yaw(i);
    for j = 1:length(sheeting_angle)
        [~,~,cT(i, j),~] = calc_objective(sheeting_angle(j));
    end
end
toc
save('cT_SA_AWA.mat', 'cT')

figure(10); clf(10); hold on;
[X,Y] = meshgrid(rad2deg(yaw), rad2deg(sheeting_angle));
surf(X, Y, cT)
c = colorbar;
c.Label.String = 'cT';

title('cT(AWA, $\delta_s$)', 'Interpreter', 'Latex')
xlabel('AWA [deg]', 'Interpreter', 'Latex');
ylabel('sheeting angle $\delta_s$ [deg]', 'Interpreter', 'Latex');