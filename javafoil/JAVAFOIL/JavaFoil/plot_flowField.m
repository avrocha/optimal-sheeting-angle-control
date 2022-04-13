function plot_flowField(awa,scale,cl,cd,cp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
awa=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('JavaFoil/flowfield.txt','r');
  line1 = fgetl(fid);
  line2 = fgetl(fid);
  line3 = fgetl(fid);
  formatSpec = '%f';
  sizeA = [5 Inf];
  data = fscanf(fid,formatSpec,sizeA)';
fclose(fid);

%data = data*scale;
x  = data(:,1);    y  = data(:,2);     vx = data(:,3);    vy = data(:,4);     v  = data(:,5);

% Reshape data
tmp = x~=x(1);
tmp = find(tmp==1);
nx  = tmp(1)-1;
xx  = x(1:nx:end);
yy  = y(1:nx);
ny  = length(x)/31;
dx  = min(diff(xx));
dy  = min(diff(yy));

figure(3);%clf;
[X,Y]  = meshgrid(min(x):dx:max(x),min(y):dy:max(y));
Vx     = reshape(vx,[nx,ny]);
Vy     = reshape(vy,[nx,ny]);
V      = reshape(v,[nx,ny]);

%contour(X,Y,V);hold on; quiver(X,Y,Vx,Vy,3);
hold on
quiver(X,Y,Vx,Vy,2.0,'b');

%cp  = [cp 0]*scale;
cp = ([cos(awa) sin(awa)])+[cp(1) 0];
%plot(cp(1),cp(2),'ro');

ct = cl*sin(awa)-cd*cos(awa);
cc = cl*cos(awa)+cd*sin(awa);
%drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'linewidth',2) ;
%drawArrow((cp(1)+[0 cd]*scale),cp(2)+[0 cl]*scale);
%drawArrow((cp(1)+[0 cd]*scale),cp(2)+[0 cl]*scale);
%plot(cp(1)+[0 cd cd 0 0]*scale,cp(2)+[0 0 cl cl 0]*scale,'k:','Linewidth',2);

end
