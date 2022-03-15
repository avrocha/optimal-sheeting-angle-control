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
global ship   counter;
diary 'Java_diary.txt'
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

ship.addRig(R1);  ship.addRig(R2);  ship.addRig(R3);  ship.addRig(R4);

scale = calc_scale();

irun=0;
for yaw = deg2rad([45 60 90 135])
    tic;
    irun=irun+1;
    figure(1);clf;figure(2);clf;figure(3);clf;
    fprintf('-------------------------------------------------------------\n');
    ship.yaw    = yaw;
    fprintf('Ship AWA = %.1f deg\n',rad2deg(ship.yaw));
    
%     if irun==1;   
      X = -ship.yaw*[1 1 1 1] + deg2rad(10);
%     end

    Mesh2Java(X);
    
    % Solve equilibrium equations using FMINCON
    ub = X+deg2rad(25);           % Upper sheet Boundary
    lb = X-deg2rad(25);           % Lower sheet Boundary
    
    calc_objective(X);
    
    counter=0;
    options  = optimoptions('fmincon','FiniteDifferenceStepSize',deg2rad(1),...
          'display','iter-detailed','OptimalityTolerance',0.001,'MaxIterations',50); % 
    [X,fval,exitflag,output] = fmincon(@(X)calc_objective(X),X,[],[],[],[],lb,ub,[],options);
    %   saveas(gca,sprintf('pix/final_geometry_awa%d',rad2deg(yaw)),'png');
    fprintf('AWA=%0.1f  Xopt=[%.1f %.1f %.1f %.1f ]\n',rad2deg(yaw), -rad2deg(X(1)),-rad2deg(X(2)),-rad2deg(X(3)),-rad2deg(X(4)));
    
    fid     = fopen([pwd,'/4_simple_wings_cl2.txt'], 'a');
    fprintf(fid,'%0.1f %.1f %.1f %.1f %.1f \n',rad2deg(yaw), rad2deg(-X(1)),rad2deg(-X(2)),rad2deg(-X(3)),rad2deg(-X(4)));
    fclose(fid);
    disp(toc);
end

diary off;
disp('done :-)');