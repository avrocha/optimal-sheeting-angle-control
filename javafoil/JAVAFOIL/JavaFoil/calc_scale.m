function scale = calc_scale();

global ship;

% Set Xval=0 to calc scale
figure(1);
X0 = zeros(1,20);
tmp = ship.yaw; ship.yaw=0;
ship.doMesh(X0,false);  % ship.doMesh(Xval,normalize,flip)
ship.yaw=tmp;
C = ship.coords;

% Handle the foil-section dividers
index = find(C(1,:)>1000); 
C(:,index)=C(:,index-1); % Put back the section dividers

xmax  = max(C(1,:));
xmin  = min(C(1,:));
ymax  = max(C(2,:));
ymin  = min(C(2,:));
ymid  = (ymin+ymax)/2;
scale = (xmax-xmin);
ship.scale=scale;
fprintf('Scale    = %.3f \n',scale);

