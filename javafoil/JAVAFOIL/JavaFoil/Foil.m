classdef Foil < handle
    % Class definition
    properties
        foilFile;
        chord;
        x,y;
        dx;
        trim;
        coords_org;
    end
    
    methods
        %-----------------------------------------------------
        function obj = Foil(foilFile,x,y,dx,chord); % Constructor method
            obj.foilFile = foilFile; % Filename : The name of the file, omittong '.txt'
            obj.chord    = chord;    % scale    : scale factor for size
            obj.x        = x;        % x        : Pivot point
            obj.y        = y;        % y        : Pivot point
            obj.dx       = dx;       % dx       : translates the scaled foil in x-dirextion
            obj.trim     = 0;        %          : rotates foil clockwise [rad]
            
            % Read the coords
           file = [pwd,'/Foils/',obj.foilFile,'.txt'];
            fid         = fopen(file);
            if fid<0; 
                disp(['Cannot open the file   ',file]);
                STOP;
            end
            streng      = fscanf(fid,'%c');
            fclose(fid);
            streng          = strrep(streng,',','.');
            obj.coords_org  = str2num(streng)';
        end
        
        %-----------------------------------------------------
        function coords = calcCoords(obj)
            coords  = [-1 0;0 1]*obj.coords_org ;    % flip profile 180deg
            coords  = coords * obj.chord;            % Scale
            coords  = coords + [obj.dx; 0 ];         % Move
            T       = [cos(obj.trim) -sin(obj.trim);sin(obj.trim) cos(obj.trim)];
            coords  = T*coords;                      % Rotate
            coords  = coords + [obj.x;obj.y];        % Move
        end
    end
end

    