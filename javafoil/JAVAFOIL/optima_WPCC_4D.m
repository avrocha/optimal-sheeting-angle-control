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
clearvars -except F*; 
clc; close all;
addpath JavaFoil; addpath Foils;  
global ship

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

% Define AWA range
% AWA = linspace(deg2rad(0), deg2rad(180), 181);
AWA = deg2rad([-120, 120, -90, 90, -60, 60, -45, 45]);
AWA = fliplr(AWA);

ship.addRig(R1);
ship.addRig(R2);
ship.addRig(R3);
ship.addRig(R4);
calc_scale();

% diary 'data\optima_calc\4D_diary.txt'

X = zeros(length(AWA), 4);
for i = 1:length(AWA)
    ship.yaw = AWA(i);

    tic;
    fprintf('4D | Ship AWA = %.1f deg\n', rad2deg(ship.yaw));
    
    if ship.yaw >= 0
        X(i, :) = -ship.yaw*[1 1 1 1] + deg2rad(10);
    else
        X(i, :) = -ship.yaw*[1 1 1 1] - deg2rad(10);
    end
    % Solve equilibrium equations using FMINCON
    ub = X(i, :) + deg2rad(25);           % Upper sheet Boundary
    lb = X(i, :) - deg2rad(25);           % Lower sheet Boundary
    
    counter=0;
    options  = optimoptions('fmincon','FiniteDifferenceStepSize',deg2rad(1),...
          'display','iter-detailed','OptimalityTolerance',0.001,'MaxIterations',50);

    [X(i, :), fval, exitflag, output] = fmincon(@(X)calc_objective(X), X(i, :), [], [], [], [], lb, ub, [], options);
    
    fprintf('AWA=%0.1f  Xopt=[%.1f, %.1f, %.1f, %.1f]\n', rad2deg(AWA(i)), rad2deg(X(i, 1)), rad2deg(X(i, 2)), rad2deg(X(i, 3)), rad2deg(X(i, 4)));
    
    toc
    
    data.X   = X;
    data.AWA = AWA;
    
%     save('data\optima_calc\4D_optima.mat', 'data')
end


% diary off;

%% Plots
load('data\optima_calc\4D_optima.mat')
fig_cnt = 1;

N = length(data.AWA);

figure(fig_cnt), clf(fig_cnt);
hold on;
plot(0:N-1, rad2deg(data.X(:, 1)))
plot(0:N-1, rad2deg(data.X(:, 2)))
plot(0:N-1, rad2deg(data.X(:, 3)))
plot(0:N-1, rad2deg(data.X(:, 4)))
xlabel('AWA [º]')
ylabel('SA [º]')
legend('Sail 1', 'Sail 2', 'Sail 3', 'Sail 4')
fig_cnt = fig_cnt + 1;