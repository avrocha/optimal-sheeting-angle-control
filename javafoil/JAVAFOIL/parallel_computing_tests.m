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

sheeting_angle_1 = deg2rad(linspace(1, 10, 10));
sheeting_angle_2 = deg2rad(linspace(11, 20, 10));
L1 = length(sheeting_angle_1);
L2 = length(sheeting_angle_2);

a = zeros(L1, L2);
b = zeros(L1, L2);
a_par = zeros(L1, L2);
b_par = zeros(L1, L2);

%% PARFOR tests
% Serial
tic
for i = 1:L1
    for j = 1:L2
        a(i,j) = fun(sheeting_angle_1(i), sheeting_angle_2(j), localShip, 1);
        b(i,j) = fun2(sheeting_angle_1(i), sheeting_angle_2(j));
    end
end
toc

% Parallel
parpool('local', 4);

tic
% ticBytes(gcp);
parfor i = 1:L1
    for j = 1:L2
        a_par(i,j) = fun(sheeting_angle_1(i), sheeting_angle_2(j), localShip, getCurrentTask().ID);
        b_par(i,j) = fun2(sheeting_angle_1(i), sheeting_angle_2(j));
    end
end
% tocBytes(gcp)
toc

poolobj = gcp('nocreate');
delete(poolobj);


%% PARFEVAL tests
p = parpool('local', 2);

% Background computation
esc_loop = @(sheeting_angle_1, sheeting_angle_2, ship, task_id) gb_loop(fun, sheeting_angle_1, sheeting_angle_2, ship, task_id);
F = parfeval(esc_loop, 1, sheeting_angle_1, sheeting_angle_2, localShip, 2);

% Main worker loop
cT = zeros(L1, L2);
for i = 1:L1
    for j = 1:L2
        cT(i,j) = fun(sheeting_angle_1(i), sheeting_angle_2(j), localShip, 1);
    end
end

cT_background = fetchOutputs(F);

poolobj = gcp('nocreate');
delete(poolobj);

%% Functions
% GB-ESC loop
function cT = gb_loop(f, sheeting_angle_1, sheeting_angle_2, ship, task_id)
    L1 = length(sheeting_angle_1);
    L2 = length(sheeting_angle_2);
    cT = zeros(L1, L2);

    for i = 1:L1
        for j = 1:L2
            cT(i,j) = f(sheeting_angle_1(i), sheeting_angle_2(j), ship, task_id);
        end
    end
end
