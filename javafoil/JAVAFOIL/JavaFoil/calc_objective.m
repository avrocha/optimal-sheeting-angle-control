function [obj,cl,CT,cd] = calc_objective(X)
    global ship counter;
    % ship.rigs(1).trim = X(1);        ship.rigs(1).foils(2).trim = -X(2);
    % ship.rigs(2).trim = X(3);        ship.rigs(2).foils(2).trim = -X(4);
    %ship.rigs(3).trim = X(5);        ship.rigs(3).foils(2).trim = -X(6);
    %ship.rigs(4).trim = X(7);        ship.rigs(4).foils(2).trim = -X(8);
    
    % Adapt to library w/ parallel computations
    process_id = '1';

    localShip = ship;
    Mesh2Java(X, localShip, process_id); 
    
    % Run JavaFoil
    Re         = 10^6;
    script_str = genJava(Re, process_id); % Generates the JavaFoil run script script.jfscript 
    javaPath   = [pwd,'/JavaFoil'];
    cmd        = ['java -cp "',javaPath,'/mhclasses.jar" -jar "', javaPath,'/javafoil.jar" Script="', script_str, '"'];
    system(cmd);

    [alfa,cl,cd,cm,cp] = readJavaResults(process_id);
     
    % Calculate lift and drag coeffs
    CT       = cl*sin(ship.yaw)-cd*cos(ship.yaw); % Propulsive force
    obj      = -CT^2*sign(CT); % Thrust
%     obj      = -cl^2;
    %fprintf('Objective function %0.4f \n',obj);
    
    % Uncomment line below to plot the flow field
%     plot_flowField(ship.yaw,ship.scale,cl,cd,cp, process_id); 
    
    hold on;zoom on
%     plot([cp cp],[-0.5 0.5],'k');%text(cp-0.02, .52,'cp','fontsize',16);
    if length(X) == 4
        title(sprintf('AWA=%0.1f  X=[%.0f/%.0f %.0f/%.0f],  cT=%.4f\n',rad2deg(ship.yaw), rad2deg(X(1)),rad2deg(X(2)),rad2deg(X(3)),rad2deg(X(4)),CT),'fontsize',16); %,rad2deg(X(5)),rad2deg(X(6)),rad2deg(X(7)),rad2deg(X(8))
    elseif length(X) == 2
        title(sprintf('AWA=%0.1f  X=[%.0f, %.0f], cT=%.4f\n',rad2deg(ship.yaw), rad2deg(X(1)), rad2deg(X(2)), CT),'fontsize',16);
    elseif length(X) == 1
        title(sprintf('AWA=%0.1f  X=%.0f,  cT=%.4f\n',rad2deg(ship.yaw), rad2deg(X), CT),'fontsize',16);
    else
       disp('Error')
       quit    
    end
   
%     drawnow;
    
    counter=counter+1;
    %saveas(gcf,['pix/bild',num2str(counter)],'png');
    fprintf('cT=%.4f, ',CT);
    fprintf('cl=%.4f, ',cl);
    fprintf('cd=%.4f, ',cd);
    fprintf('cm=%.4f, ',cm);
    fprintf('cp=%.4f, ',cp);
    if length(X) == 4
        fprintf('X=[%.0f %.0f %.0f %.0f ]\n',rad2deg(X(1)),rad2deg(X(2)),rad2deg(X(3)),rad2deg(X(4)));
    elseif length(X) == 2
        fprintf('X=[%.0f %.0f]\n', rad2deg(X(1)), rad2deg(X(2)));
    elseif length(X) == 1
        fprintf('X=%.0f\n', rad2deg(X));
    else
       disp('Error')
       quit    
    end
end



