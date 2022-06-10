
function out = calc_objective_mod(X, ship, process_id)
    % Inputs:
    % X          : [1 x n] sheeting angle vector; X(n) corresponds to the foremost wingsail sheeting angle
    % ship       : [1 x 1] Ship object.
    % process_id : [1 x 1] (int) ID of current task (process) to avoid shared memory in multiprocessing

    if nargin < 3
        process_id = '1';
    else    
        process_id = sprintf('%d', int8(process_id));
    end

    % Generates Java-scale mesh-file
    Mesh2Java(X, ship, process_id); 
    
    % Run JavaFoil
    Re         = 10^6;
    script_str = genJava(Re, process_id); % Generates the JavaFoil run script script.jfscript 
    javaPath   = [pwd,'/JavaFoil'];
    cmd        = ['java -cp "',javaPath,'/mhclasses.jar" -jar "', javaPath,'/javafoil.jar" Script="', script_str, '"'];
    system(cmd);
    
    [~, cL, cD, ~, cP] = readJavaResults();
     
    % Calculate lift and drag coeffs
    cT       = cL * sin(ship.yaw) - cD * cos(ship.yaw); % Propulsive force
    out.obj  = -cT^2 * sign(cT); % Thrust
    % out.obj  = -cl^2;
    out.cT  = cT;
    out.cL  = cL;
    out.cD  = cD;

    % Uncomment line below to plot the flow field
%     plot_flowField(ship.yaw,ship.scale,cL,cD,cP); 
end
