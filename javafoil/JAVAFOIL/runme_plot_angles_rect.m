clear; %clc; close all;  
addpath JavaFoil; addpath Foils;  
global ship   counter;


ship   = Ship(200); % (yaw,Pivot )Create the ship object
Cw     = 25; % Wing chord
Cf     = 12.5; % Flap chord
R1     = Rig(26,-35/2); % pivot x,y,  
R1.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord,
% R1.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

R2     = Rig(26,35/2); % pivot x,y, 
R2.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% R2.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord,

R3     = Rig(26+150,-35/2); % pivot x,y,     
R3.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% R3.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

R4     = Rig(26+150,35/2); % pivot x,y,  
R4.addFoil(Foil('NACA0018',0,0,6.25,Cw)); % foilFile,x,y,dx,chord
% R4.addFoil(Foil('NACA0009',-5,0,0,Cf)); % foilFile,x,y,dx,chord

ship.addRig(R1);  ship.addRig(R2);  ship.addRig(R3);  ship.addRig(R4);

scale = calc_scale();

anglesWings = csvread('trim_angles_4wings_cl_rectangle2_conf.txt');


yawVec =deg2rad(anglesWings(:,1));
for ii = 1:1:length(yawVec)
    yaw=yawVec(ii);
    
    ship.yaw    = yaw;
    fprintf('Ship AWA = %.1f deg\n',rad2deg(ship.yaw));
    
    X=-deg2rad(anglesWings(ii,2:end));
%     X(1)=X(1)+deg2rad(20);
%     X(3)=X(3)+deg2rad(20);
    
%     ship.doMesh(X,true); 
%     saveas(gca,sprintf('pix/rectangle_awa%i',rad2deg(yaw)),'png');
    
    figure(ii);
    Mesh2Java(X);
    calc_objective(X);
    saveas(gca,sprintf('pix/flow_rectangle_awa%i',(anglesWings(ii,1))),'png');
%     ship.doMesh(X,true);
end
