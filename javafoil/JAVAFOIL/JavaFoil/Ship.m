classdef Ship < handle
    % Class definition
    properties
        rigs;
        yaw ;
        length;
        x,y;  % Pivot point
        myFrame=0;
        myFilm ;
        scale
        x_mid;
        y_mid;
        coords=[]; % Preliminary set of coords to export
    end
    
    methods
        %-----------------------------------------------------
        function obj = Ship(length); % Constructor method
            obj.yaw    = 0;
            obj.length = length;
            obj.x      = 0.0;
            obj.y      = 0.0;
        end
        %-----------------------------------------------------
        function addRig(obj,rig); 
           obj.rigs = [obj.rigs rig];
        end
        %-----------------------------------------------------
        function coords = transform(obj,coords)
            T = [cos(obj.yaw) -sin(obj.yaw);sin(obj.yaw) cos(obj.yaw)];
            coords = T*coords + [obj.x;obj.y];
        end
        %-----------------------------------------------------
        function doMesh(obj,Xval, hull); 
          C=[];
          x_min = 9999999;x_max = -9999999;y_min = 9999999;y_max = -9999999;
          Tawa  = [cos(obj.yaw) -sin(obj.yaw);sin(obj.yaw) cos(obj.yaw)];
          ixval = 0;
          for ir=1:length(obj.rigs);
              rig      = obj.rigs(ir); 
              rig.trim = Xval(ir);%2*ir-1);
              for ifoil=1:length(rig.foils);
                foil     = rig.foils(ifoil);
%                 if (mod(ifoil,2)==0); foil.trim = Xval(2*(ir-1)+ifoil);else foil.trim=0;end %
                coords = foil.calcCoords;        % 1) Rotate with foil trim
                coords = rig.transform(coords);  % 2) Rotate with rig trim
                coords = obj.transform(coords);  % 3) Rotates with ship awa
                %rig.foils(ifoil).coords=coords;  % save  for later
              
                C = [C ,coords,[999999.9;999999.9]]; % To show Javafoil end of profile
               end
          end
          obj.coords=C;
%           obj.plotCoords;
          hold on;
          if hull; obj.plotHull;end
          hold on;
          
        end
        %-----------------------------------------------------
        function plotCoords(obj)
           C =  obj.coords(:,1:end-1);
           C(:,C(1,:)>1000)=0*C(:,C(1,:)>1000); % Fix the 999999
           plot(C(1,:)./obj.scale,C(2,:)./obj.scale,'.');axis equal;grid on;
        end
        %-----------------------------------------------------
        function saveToJavaFoil(obj)
           C =  obj.coords(:,1:end-1); % Take away the last 9999 of the last foil
           fid     = fopen([pwd,'/JavaFoil/Mesh_fullsize.txt'], 'w');
           fprintf(fid, '%12.6f %12.6f\n', C);
           fclose(fid);
        end
        %-----------------------------------------------------
        function plotHull(obj); 
           T  = [cos(obj.yaw) sin(obj.yaw);-sin(obj.yaw) cos(obj.yaw)];
%            boat=(1.2*[ -0.2 0;0.1 0.1;1 0.1;1 -0.1;0.1 -0.1;-0.2 0]-[0.5 0])*Tawa+[0.5 0];
%            boat=(1.2*T*([ -0.05 -0.1;-0.05 0.1;1 0.1; 1.1 0;1 -0.1;-0.05 -0.1]-[0.5 0.5;0.5 0.5;0.5 0.5;0.5 0.5;0.5 0.5;0.5 0.5])')-[0. 0.5]';%-[0.5 0]');%+[0.5 -0.0]';
           boat=-1.5*T*([ -0.05 -0.1;-0.05 0.1;1 0.1; 1.1 0;1 -0.1;-0.05 -0.1]'+[-0.5 0]')+[0.5 0]';
           boat = boat*obj.length/obj.scale;
           p1= plot(boat(1,:),boat(2,:),':k','LineWidth',2);
           p1.Color(4) = 0.7;
           %text(boat(1,1)-0.1,boat(2,1)+0.3,['AWA=',num2str(rad2deg(obj.yaw),'%5.1f'),'^o'],'fontsize',12);
           title(['AWA=',num2str(rad2deg(obj.yaw),'%5.1f'),'^o']);
           axis equal;
           grid on;zoom on;
%            set(gca,'Ydir','reverse');
           %xlabel('x');ylabel('y');
        end
        %-----------------------------------------------------
        function analyseResults(obj)
           
        end
    end
end

    