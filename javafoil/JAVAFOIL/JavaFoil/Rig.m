classdef Rig < handle
    % W Rig consists of one or more foils
    properties
        x, y ;   % Pivot point = mast
        trim ;
        foils;  % The wing and the flap (treated the same)
    end
    
    methods
        %-----------------------------------------------------
        function obj = addFoil(obj,foil); 
           obj.foils = [obj.foils foil];
        end
        %-----------------------------------------------------
        function obj = Rig(x,y); % Constructor method
            obj.x      = x;
            obj.y      = y;
            obj.trim   = 0;
        end
        %-----------------------------------------------------
        function coords = transform(obj,coords) 
            T      = [cos(obj.trim) -sin(obj.trim);sin(obj.trim) cos(obj.trim)];
            coords = T*coords+[obj.x;obj.y];         
        end
        %-----------------------------------------------------
    end  
end
    