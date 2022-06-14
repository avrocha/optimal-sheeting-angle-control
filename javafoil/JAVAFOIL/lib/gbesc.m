function [u, y, dy, y_hat] = gbesc(ship, J, dt, N, f, A, fc_hp, fc_lp, K, u0, AWA, lp_bool, cT_filter, cT_filter_param, FF)
    % Inputs:
    % - J        : optimization criterion [function handle]
    % - dt       : simulation step [s]
    % - N        : simulation lenght 
    % - f        : sinusoidal dithers frequency [Hz]
    % - A        : sinusoidal dithers amplitude [rad]
    % - fc_hp    : HPF cut-off frequency [Hz]
    % - fc_lp    : LPPF cut-off frequency [Hz]
    % - K        : integrator gain
    % - u0       : initial sheeting angle [rad]
    % - AWA      : [1 x N] time variant AWA
    % - lp_bool  : use LPF [boolean]
    % - cT_filter: cT filter type {'RAW', 'EMA', 'LPF'}
    % - cT_filter_param : Parameter for cT filter (0, alpha, cut-off frequency, respectively)
    % - FF       : [1 x N] time variant FF
    
    % Outputs:
    % - u : control variable
    % - y : criterion output
    % - dy: criterion gradient estimate
    
    if (length(A) ~= length(f))
        fprintf('Dimensions do not match\n.')
        return       
    end    
    
    % Data structures
    n     = length(f);
    u_hat = zeros(n, N+1);
    u     = [u0, zeros(n, N)];
    y     = zeros(1, N);
    y_hat = zeros(1, N);
    dy    = zeros(n, N);
    hpf   = zeros(1, N);
    zeta  = zeros(n, N);
    lpf   = zeros(n, N);

    % cT filtering
    switch cT_filter
        case 'RAW'
            disp('cT filtering: No cT filter selected.')

        case 'EMA'
            disp('cT filtering: EMA cT filter selected (select rate parameter frequency manualy).')
            alpha = cT_filter_param;

        case 'LPF'
            disp('cT filtering: LPF cT filter selected (select cut-off frequency manualy).')
            bworder        = 5;
            [bl_y ,al_y]   = butter(bworder, cT_filter_param*dt, 'low');
            bl_y = fliplr(bl_y);
            al_y = fliplr(al_y);
            M_y  = bworder + 1;

        otherwise
            disp('cT filtering: Invalid filter selection.')
    end
    
    % Filters
    bworder = 5;
    % HPF
    [bh,ah]   = butter(bworder, fc_hp*dt, 'high');
    bh = fliplr(bh);
    ah = fliplr(ah);

    % LPF
    [bl,al]   = butter(bworder, fc_lp*dt, 'low');
    bl = fliplr(bl);
    al = fliplr(al);

    M = bworder + 1; % filter length

    
    if (length(A) ~= length(f)) && (length(f) ~= n) && (n ~= length(ship.yaw()))
        fprintf('Dimensions do not match\n.')
        return       
    end

    tic
    for i = 1:N
        if rem(i, 10) == 0
            fprintf("Iteration %d\n", i);
        end
        t = i*dt;
   
        ship.yaw = AWA(i);
        
        y(i) = J(u(:, i), ship);
        
        % cT smoothing
        switch cT_filter
            case 'RAW'
                y_hat(i) = y(i);

            case 'LPF'
                if i > M_y
                    y_hat(i) = 1/al_y(end) * (-al_y(1:end-1) * y_hat(i-M_y+1:i-1)' + bl_y * y(i-M_y+1:i)');
                else
                    y_hat(i) = y(1); % Start by filling filter
                end

            case 'EMA'
                if i == 1
                    y_hat(i) = y(i);
                else
                    y_hat(i) = (1-alpha)*y_hat(i-1) + alpha*y(i);
                end

        end 

        %HPF
        if i >= M
            hpf(i) = 1/ah(end) * (-ah(1:end-1) * hpf(i-M+1:i-1)' + bh * y_hat(i-M+1:i)');
        else
            for j = 1:i
                hpf(i) = hpf(i) + bh(end-j+1)*y_hat(i-j+1);
            end
    
            for j = 2:i
                hpf(i) = hpf(i) - ah(end-j+1)*hpf(i-j+1);
            end
            
            hpf(i) = 1/ah(end) * hpf(i);
        end
        
        if lp_bool
            zeta(:, i) = hpf(i) * sin(2*pi*f*t);
            
            % LPF
            if i >= M              
                lpf(:, i) = 1/al(end) * (-al(1:end-1) * lpf(:, i-M+1:i-1)' + bl * zeta(:, i-M+1:i)');
            else
                for j = 1:i
                    lpf(:, i) = lpf(:, i) + bl(end-j+1) .* zeta(:, i-j+1);
                end
        
                for j = 2:i
                    lpf(:, i) = lpf(:, i) - al(end-j+1) .* lpf(:, i-j+1);
                end
                
                lpf(:, i) = 1/al(end) .* lpf(:, i);
            end
    
            dy(:, i) = lpf(:, i);
        
        else
            dy(:, i) = hpf(i) * sin(2*pi*f*t);

        end
        
        % Parameter estimate - Single Integrator
        u_hat(:, i+1) = u_hat(:, i) + dt * K * dy(:, i); % single integrator
        
        % Add dither
        u(:, i+1)     = FF(i) + u_hat(:, i+1) + A .* sin(2*pi*f*t);
        
        % Error condition
        if any(u(i+1) > pi) || any(u(i+1) < -pi) 
            break
        end
    end
    toc
end