function out = calc_objective_mod(X)
    global ship counter;
    
    Mesh2Java(X); % Generates Java-scale mesh-file
    
    awa = 0;  % For JavaFoil
    Re  = 10^6;
    genJava(Re, awa); % Generates the JavaFoil run script script.jfscript 
    javaPath = [pwd,'/JavaFoil'];
    cmd = ['java -cp "',javaPath,'/mhclasses.jar" -jar "', javaPath,'/javafoil.jar" Script="',javaPath,'/script.jfscript"'];
    system(cmd);     % Excecute the JavaFoil calculations
    
    [~, cL, cD, ~, ~] = readJavaResults();
     
    % Calculate lift and drag coeffs
    cT       = cL * sin(ship.yaw) - cD * cos(ship.yaw); % Propulsive force
    out.obj  = -cT^2 * sign(cT); % Thrust
    % out.obj  = -cl^2;
    out.cT  = cT;
    out.cL  = cL;
    out.cD  = cD;
end
