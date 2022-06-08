function script_str = genJava(Re, process_id)
    %======================================================================
    % Generates input file for JavaFoil
    %======================================================================
    % Inputs:
    % Re         : [1 x 1] Reynolds number
    % process_id : [1 x 1] (str) ID of current task (process) - suffix for W/R files
        
    if ~ischar(process_id)
        disp('genJava: argument process_id is not valid.\n');
        return;
    end
    
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
    script_str = strcat('JavaFoil/script_', process_id, '.jfscript');
    fid        = fopen(script_str, 'w');
    
    % Header general information
    line{1} = (['// recorded on ',datestr(now),' by ',getenv('USER')]);
    line{2} = ['Options.Country(',Country,')']; % Settings
    line{3} = ['Options.TransitionModel(',TransitionModel,')']; % Settings
    line{4} = ['Options.StallModel(',StallModel,')']; % Settings
    line{5} = ['Geometry.Open("',pwd,'/JavaFoil/Mesh_',process_id,'.txt")']; % Generate foil geometry
    line{6} = ['Polar.Analyze(',ReFirst,';',ReLast,';',ReDelta,';',AngleFirst,';',AngleLast,';',AngleDelta,';',TransUpper,';',TransLower,';',Roughness,';',AddToPlotsFlag,')']; % Analyze foil
    
    if (ismac);line{6} = strrep(line{6},';',':');end 
   
    % Output files defined per each processor
    line{7} = ['Polar.Save(',pwd,'/JavaFoil/output_',process_id,'.txt)'];
    line{8} = ['FlowField.Save(',pwd,'/JavaFoil/flowfield_',process_id,'.txt)'];
    line{9} = ['Exit()'];
    % Write file
    for i=1:length(line)
        fprintf(fid,'%s\n',line{i});
    end
    
    % Close file return to main
    fclose(fid);

end
