function cT = interp_criterion_irregular(xAWA, xSA, V, x, interp_type, f, ship)
    % Interpolate irregular grid: 2-step interpolation (1st interpolate in delta, for AWA
    % neighbors; 2nd interpolate in AWA)
    % --- 
    % xAWA        : AWA axis
    % xSA         : SA irregular axes
    % V           : cT data
    % x           : Query point [AWA, Sheeting angles]
    % interp_type : Interpolation type: 'linear' | 'poly'
    % f           : Function to call when sheeting angle is outside boundaries
    % --
    % Outputs
    % cT          : Interpolated thrust coefficient at (x, AWA)
    
    if ~strcmp(interp_type, 'linear') && ~strcmp(interp_type, 'cubic')
        disp('interp_type is not valid. \n')
    end

    n = size(xSA, 1);
    if n ~= 4
        disp('Criterion is not developed for configurations different 4 wingsails.\n')
    end
    
%     fprintf("Inputs: AWA = %f, SA = %s\n", rad2deg(x(1)), num2str(rad2deg(x(2:end))));

    % Get AWA neighbor levels 
    [~, iAWA] = min(abs(x(1) - xAWA));
    
    if x(1) >= xAWA(iAWA) % x(1) near left neighbour
        iiAWA = max(1, iAWA):min(length(xAWA), iAWA+1);
    else % x(1) near right neighbour
        iiAWA = max(1, iAWA-1):min(length(xAWA), iAWA);
    end

%     fprintf("AWA neighbors: %f | %f\n", rad2deg(xAWA(iiAWA(1))), rad2deg(xAWA(iiAWA(2))))

    out_of_bounds = 0;
    for i = 1:n
        % Check boundary conditions on AWA neighbors

        % Max min of each neighbor
        lb = max([min(xSA(i, :, iiAWA(1)), [],'all'), min(xSA(i, :, iiAWA(2)), [],'all')]);
        % Min max of each neighbor
        ub = min([max(xSA(i, :, iiAWA(1)), [],'all'), max(xSA(i, :, iiAWA(2)), [],'all')]);
        
        out_of_bounds = out_of_bounds | (x(i+1) < lb) | (x(i+1) > ub);
        
        if out_of_bounds
            fprintf("x(%d) %f | lb = %f | ub = %f\n", i, rad2deg(x(i+1)), rad2deg(lb), rad2deg(ub));
            fprintf("OUT = %d\n", out_of_bounds)
            disp('\tInterpolated criterion is out of bounds\n.')
            break; 
        end
    end

    
    X = cell(5, 1);    
    if ~out_of_bounds
        % Inf neighbor
        [X{1}, X{2}, X{3}, X{4}] = ndgrid(squeeze(xSA(1, :, iiAWA(1))), ...
                                            squeeze(xSA(2, :, iiAWA(1))), ...
                                            squeeze(xSA(3, :, iiAWA(1))), ...
                                            squeeze(xSA(4, :, iiAWA(1))));
        

        cT_inf = interpn(X{1}, X{2}, X{3}, X{4}, squeeze(V(iiAWA(1), :, :, :, :)), x(2), x(3), x(4), x(5), interp_type);

        % Sup neighbor
        [X{1}, X{2}, X{3}, X{4}] = ndgrid(squeeze(xSA(1, :, iiAWA(2))), ...
                                            squeeze(xSA(2, :, iiAWA(2))), ...
                                            squeeze(xSA(3, :, iiAWA(2))), ...
                                            squeeze(xSA(4, :, iiAWA(2))));

    
        cT_sup = interpn(X{1}, X{2}, X{3}, X{4}, squeeze(V(iiAWA(2), :, :, :, :)), x(2), x(3), x(4), x(5), interp_type);
        
        cT = interp1(xAWA(iiAWA), [cT_inf, cT_sup], x(1));

%         fprintf("cT_inf = %f | cT_sup = %f | cT_interp = %f\n--\n", cT_inf, cT_sup, cT);

    else
        cT = f(x(2:end), ship);
%         fprintf("cT %f\n--\n", cT);

    end

end