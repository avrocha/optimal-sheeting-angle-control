function [alfa,cl,cd,cm,cp] = readJavaResults()

% Open results file from JavaFoil
filename = [pwd,'/JavaFoil/output.txt'];
fid = fopen(filename,'r');
Intro = textscan(fid,'%s',5,'Delimiter','\n');

data = fscanf(fid,'%f',[11 Inf])'; %  Read the data
fclose(fid);

% Extract cl and cd and return to main
alfa = deg2rad(data(:,1)); % Förihelvete
cl   = data(:,2);
cd   = data(:,3);
cm   = data(:,4);
cp   = data(:,11);

