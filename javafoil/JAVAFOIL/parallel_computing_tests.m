addpath JavaFoil; addpath Foils; addpath lib
global ship;

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

AWA = deg2rad(45);
localShip = ship;
localShip.yaw = AWA;

fun = @(delta_1, delta_2, localShip, task_id) getfield(calc_objective_mod([delta_1, delta_2], ...
                                                localShip, task_id), 'cT');
fun2 = @(delta_1, delta_2) delta_1 * delta_2;

sangle_1 = deg2rad(linspace(1, 10, 10));
sangle_2 = deg2rad(linspace(11, 20, 10));
Ldelta_1 = length(sangle_1);
Ldelta_2 = length(sangle_2);

a = zeros(Ldelta_1, Ldelta_2);
b = zeros(Ldelta_1, Ldelta_2);

a_par = zeros(Ldelta_1, Ldelta_2);
b_par = zeros(Ldelta_1, Ldelta_2);

%% Serial
tic
for i = 1:Ldelta_1
    for j = 1:Ldelta_2
        a(i,j) = fun(sangle_1(i), sangle_2(j), localShip, 1);
        b(i,j) = fun2(sangle_1(i), sangle_2(j));
    end
end
toc

%% Parallel
parpool('local', 4);

tic
% ticBytes(gcp);
parfor i = 1:Ldelta_1
    d1 = sangle_1(i);
    for j = 1:Ldelta_2
        a_par(i,j) = fun(d1, sangle_2(j), localShip, getCurrentTask().ID);
        b_par(i,j) = fun2(d1, sangle_2(j));
    end
end
% tocBytes(gcp)
toc

poolobj = gcp('nocreate');
delete(poolobj);

