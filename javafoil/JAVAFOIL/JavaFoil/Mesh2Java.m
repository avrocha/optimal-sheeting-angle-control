function [] = Mesh2Java(Xval);
global ship;

% Now mesh with actual Xval
% figure(2);
ship.doMesh(Xval,false);  % ship.doMesh(Xval,normalize,flip)
C = ship.coords;
index = find(C(1,:)>1000); 
for ii=1:length(index);C(:,index(ii))= C(:,1);end
xmax  = max(C(1,:));
xmin  = min(C(1,:));
ymax  = max(C(2,:));
ymin  = min(C(2,:));
ymid  = (ymin+ymax)/2;

C(1,:) = (C(1,:)-xmin)/ship.scale;
C(2,:) = (C(2,:)-ymid)/ship.scale;
C(1,:) = -C(1,:)+1;

% Uncomment lines below to plot the foils
figure(3);clf;
plot(C(1,:),C(2,:),'.');
grid on;axis equal;hold on

C(:,index)=ones(2,length(index))*99999999; % Put back the section dividers

% Save the mesh for JavaFoil
fid     = fopen('JavaFoil/Mesh.txt', 'w');
fprintf(fid, '%12.6f %12.6f\n', C(:,1:end-1));
fclose(fid);

