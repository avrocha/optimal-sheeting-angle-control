function genJava(Re, awa)
%==========================================================================
% Generates input file for JavaFoil
%==========================================================================

% Set options
Country          = '0';                 % Units and such
TransitionModel  = '1';                 % Extended Eppler see reference manual
StallModel       = '0';         

% Set analysis conditions
ReFirst    = num2str(round(Re));
ReLast     = num2str(round(Re));
ReDelta    = num2str(round(Re));                % Reynolds number
AngleFirst = num2str(0);  % 20 degrees
AngleLast  = num2str(0); % 80 degrees
AngleDelta = '1';       %num2str(Alpha(2)-Alpha(1)); % Step in alpha
TransUpper = '100';               % Default
TransLower = '100';               % Default
Roughness  = '1';                 % Surface roughness, default=smooth=0
AddToPlotsFlag = '0';             % Default

% Create file using filepath
fid = fopen('JavaFoil/script.jfscript','w');

% Header general information
line{1} = (['// recorded on ',datestr(now),' by ',getenv('USER')]);
line{2} = ['Options.Country(',Country,')']; % Settings
line{3} = ['Options.TransitionModel(',TransitionModel,')']; % Settings
line{4} = ['Options.StallModel(',StallModel,')']; % Settings
line{5} = ['Geometry.Open("',pwd,'/JavaFoil/Mesh.txt")']; % Generate foil geometry
line{6} = ['Polar.Analyze(',ReFirst,';',ReLast,';',ReDelta,';',AngleFirst,';',AngleLast,';',AngleDelta,';',TransUpper,';',TransLower,';',Roughness,';',AddToPlotsFlag,')']; % Analyze foil
if (ismac);line{6} = strrep(line{6},';',':');end 
line{7} = ['Polar.Save(',pwd,'/JavaFoil/output.txt)'];
line{8} = ['FlowField.Save(',pwd,'/JavaFoil/flowfield.txt)'];
line{9} = ['Exit()'];
% Write file
for i=1:length(line)
    fprintf(fid,'%s\n',line{i});
end

% Close file return to main
fclose(fid);
