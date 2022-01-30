function [obj cl CT] = calc_objective(X)

global ship counter;
%ship.rigs(1).trim = X(1);        ship.rigs(1).foils(2).trim = -X(2);
%ship.rigs(2).trim = X(3);        ship.rigs(2).foils(2).trim = -X(4);
%ship.rigs(3).trim = X(5);        ship.rigs(3).foils(2).trim = -X(6);
%ship.rigs(4).trim = X(7);        ship.rigs(4).foils(2).trim = -X(8);

Mesh2Java(X); % Generates Java-scale mesh-file

awa = 0;  % For JavaFoil
Re  = 10^6;
genJava(Re, awa); % Generates the JavaFoil run script script.jfscript 
javaPath = [pwd,'/JavaFoil'];
cmd = ['java -cp "',javaPath,'/mhclasses.jar" -jar "', javaPath,'/javafoil.jar" Script="',javaPath,'/script.jfscript"'];
system(cmd);     % Excecute the JavaFoil calculations
[alfa,cl,cd,cm,cp] = readJavaResults();
 
% Calculate lift and drag coeffs
CT       = cl*sin(ship.yaw)-cd*cos(ship.yaw); % Propulsive force
% obj      = -CT^2*sign(CT); % Thrust
obj      = -cl^2;
%fprintf('Objective function %0.4f \n',obj);
plot_flowField(ship.yaw,ship.scale,cl,cd,cp); 

% hold on;zoom on
% plot([cp cp],[-0.5 0.5],'k');%text(cp-0.02, .52,'cp','fontsize',16);
title(sprintf('AWA=%0.1f  X=[%.0f/%.0f %.0f/%.0f],  cT=%.4f\n',rad2deg(ship.yaw), rad2deg(X(1)),rad2deg(X(2)),rad2deg(X(3)),rad2deg(X(4)),CT),'fontsize',16); %,rad2deg(X(5)),rad2deg(X(6)),rad2deg(X(7)),rad2deg(X(8))
drawnow;

counter=counter+1;
%saveas(gcf,['pix/bild',num2str(counter)],'png');

fprintf('cT  = %.4f,  ',CT);
fprintf('cl=%.4f, ',cl);
fprintf('cd=%.4f, ',cd);
fprintf('cm=%.4f, ',cm);
fprintf('cp=%.4f, ',cp);
fprintf(' X=[%.0f %.0f %.0f %.0f ]\n',rad2deg(X(1)),rad2deg(X(2)),rad2deg(X(3)),rad2deg(X(4)));%,rad2deg(X(5)),rad2deg(X(6)),rad2deg(X(7)),rad2deg(X(8)));
%%.0f %.0f %.0f %.0f 

end



