function [u, y, dy] = gbesc(ship, J, dt, N, f, A, fc_hp, fc_lp, K, u0, AWA, lp_bool)
    % Inputs:
    % - J      : optimization criterion [function handle]
    % - dt     : simulation step [s]
    % - N      : simulation lenght 
    % - f      : sinusoidal dithers frequency [Hz]
    % - A      : sinusoidal dithers amplitude [rad]
    % - fc_hp  : HPF cut-off frequency [Hz]
    % - fc_lp  : LPPF cut-off frequency [Hz]
    % - K      : integrator gain
    % - u0     : initial sheeting angle [rad]
    % - AWA    : time variant AWA
    % - lp_bool: use LPF [boolean] 
    % Outputs:
    % - u : control variable
    % - y : criterion output
    % - dy: criterion gradient estimate
    
    bworder = 5;
    % HPF
    [bh,ah]   = butter(bworder, fc_hp*dt, 'high');
    % LPF
    [bl,al]   = butter(bworder, fc_lp*dt, 'low');

    % Data structures
    n     = length(f);
    u_hat = [u0, zeros(n, N)];
    u     = [u0, zeros(n, N)];
    y     = zeros(1, N);
    dy    = zeros(n, N);
    hpf   = zeros(1, N);
    zeta  = zeros(n, N);
    lpf   = zeros(n, N);
    
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
        % Avoid numerical singularities
        if i > 1 && (y(i) > 1.5*y(i-1) || y(i) < 0.5*y(i-1))
                y(i) = y(i-1);
        end
    
        if i >= bworder+1
            for j = 1:bworder+1
                hpf(i) = hpf(i) + bh(j)*y(i-j+1);
            end
    
            for j = 2:bworder+1
                hpf(i) = hpf(i) - ah(j)*hpf(i-j+1);
            end
            
            hpf(i) = 1/ah(1) * hpf(i);
        else
            for j = 1:i
                hpf(i) = hpf(i) + bh(j)*y(i-j+1);
            end
    
            for j = 2:i
                hpf(i) = hpf(i) - ah(j)*hpf(i-j+1);
            end
            
            hpf(i) = 1/ah(1) * hpf(i);
        end
        
        if lp_bool
            zeta(:, i) = hpf(i) * sin(2*pi*f*t);
            
            % LPF
            if i >= bworder+1
                for j = 1:bworder+1
                    lpf(:, i) = lpf(:, i) + bl(j) .* zeta(:, i-j+1);
                end
        
                for j = 2:bworder+1
                    lpf(:, i) = lpf(:, i) - al(j) .* lpf(:, i-j+1);
                end
                
                lpf(:, i) = 1/al(1) .* lpf(:, i);
            else
                for j = 1:i
                    lpf(:, i) = lpf(:, i) + bl(j) .* zeta(:, i-j+1);
                end
        
                for j = 2:i
                    lpf(:, i) = lpf(:, i) - al(j) .* lpf(:, i-j+1);
                end
                
                lpf(:, i) = 1/al(1) .* lpf(:, i);
            end
    
            dy(:, i) = lpf(:, i);
        
        else
            dy(:, i) = hpf(i) * sin(2*pi*f*t);

        end
        
        % Parameter estimate - Single Integrator
        u_hat(:, i+1) = u_hat(:, i) + dt * K * dy(:, i); % single integrator
        
        % Add dither
        u(:, i+1)     = u_hat(:, i+1) + A .* sin(2*pi*f*t);
        
        % Error condition
        if any(u(i+1) > pi) || any(u(i+1) < -pi) 
            break
        end
    end
    toc
end