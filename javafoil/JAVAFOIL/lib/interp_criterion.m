function cT = interp_criterion(X, V, x, interp_type, f, ship)
    % Inputs:
    % X           : [n x 1] Cell array of data grids
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

    n = size(X, 1);
    
    out_of_bounds = 0;
    for i = 2:n
        out_of_bounds = out_of_bounds | (x(i) < min(X{i}, [],'all')) | (x(i) > max(X{i}, [],'all'));
        
        if out_of_bounds
            disp('\tInterpolated criterion is out of bounds\n.')
            break; 
        end
    end
    
    if ~out_of_bounds
        switch n
            case 2 % AWA + 1 Wing
                cT = interpn(X{1}, X{2}, V, x(1), x(2), interp_type);
            case 3 % AWA + 2 Wings
                cT = interpn(X{1}, X{2}, x{3}, V, x(1), x(2), x(3), interp_type);
            case 5 % AWA + 4 Wings
                cT = interpn(X{1}, X{2}, X{3}, X{4}, X{5}, V, x(1), x(2), x(3), x(4), x(5), interp_type);
            otherwise
                disp('Criterion is not developed for configurations different than 1/2/4 wingsails.\n')
                return
        end
    else
        cT = f(x(2:end), ship);
    end

end