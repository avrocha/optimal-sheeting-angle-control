function [alfa,cl,cd,cm,cp] = JavaCalc(Foil,Re,awa)

global javaPath ;
%==========================================================================
% JavaCalc
% 
% FoilCalc is a tool for calculating lift and drag coefficients for 2-d
% general 2-d foils. The script uses JavaFoil and detailed documentation 
% of JavaFoil may be found here: http://www.mh-aerotools.de/airfoils/javafoil.htm
%
%==========================================================================


% Generate JavaFoil inputfile
genJava(Re, Foil, awa);

% Run analysis
system(['java -cp "',javaPath,'mhclasses.jar" -jar "',javaPath,'javafoil.jar" Script="',javaPath,'script.jfscript"']);

% Read results
[alfa,cl,cd,cm,cp] = readJava();


