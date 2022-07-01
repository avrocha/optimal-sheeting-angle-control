function yfilt= medfilt4(y, windowSize)
% Median filter for 4D structure
% Boundary points are not filtered for simplicity of the algorithm
%   Inputs
%   y                     : 4D structure
%   windowSize [optional] : 4D vector with window size per dimension
%   --
%   Outputs
%   yfilt   : 4D filtered structure
    
    % Error conditions
    if length(size(y)) ~= 4
        disp('Input "y" is not 4D.')
            return
    end

    if nargin > 2
        if length(windowSize) ~= 4
            disp('Argument "windowSize" should be of length 4.')
            return
        end

        if any(rem(windowSize, 2) == 0)
            disp('Argument "windowSize" invalid: Window size must be an odd natural number.')
            return
        end
        
        % Make windowSize row oriented
        if size(windowSize, 1) > size(windowSize, 2)
            windowSize = windowSize';
        end

    else
        windowSize = 3 * ones(4, 1);

    end

    if any(size(y) < windowSize)
        disp('Argument "windowSize" invalid: Window size must be smaller than data dimension.')
            return
    end
    
    % Structure init
    yfilt      = y;
    windowEdge = floor(windowSize ./ 2);
    
    % Loop through interior points
    for iw = 1 + windowEdge(1):size(y, 1) - windowEdge(1)
        for ix = 1 + windowEdge(2):size(y, 2) - windowEdge(2)
            for iy = 1 + windowEdge(3):size(y, 3) - windowEdge(3)
                for iz = 1 + windowEdge(4):size(y, 4) - windowEdge(4)                    
                    % Filtering window
                    nDwindow = y(iw-windowEdge(1):iw+windowEdge(1),...
                                  ix-windowEdge(2):ix+windowEdge(2),...
                                   iy-windowEdge(3):iy+windowEdge(3),...
                                    iz-windowEdge(4):iz+windowEdge(4));
                    
                    % Compute median (ignore NaN)
                    yfilt(iw, ix, iy, iz) = median(nDwindow, 'all', 'omitnan');                        

                end
            end
        end
    end   

end