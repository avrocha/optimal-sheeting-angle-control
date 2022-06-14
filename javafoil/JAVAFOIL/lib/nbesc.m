function [u, y, dy, ddy, ddy_inv, y_hat] = nbesc(ship, J, dt, N, f, A, fc_hp, fc_lp, K, u0, wric, ddy0, AWA, lp_bool, cT_filter, cT_filter_param, FF)
    % Inputs:
    % - J         : optimization criterion [function handle]
    % - dt        : simulation step [s]
    % - N         : simulation lenght 
    % - f         : sinusoidal dither frequency [Hz]
    % - A         : sinusoidal dither amplitude [rad]
    % - fc_hp     : HPF cut-off frequency [Hz]
    % - fc_lp     : LPF cut-off frequency [Hz]
    % - K         : integrator gain
    % - u0        : initial sheeting angle [rad]
    % - wric      : Ricatti filter parameter
    % - ddy0      : hessian inverse initial value
    % - AWA       : time variant AWA
    % - lp_bool   : use LPF [boolean]
    % - cT_filter : cT filter type {'RAW', 'EMA', 'LPF'}
    % - cT_filter_param : Parameter for cT filter (0, alpha, cut-off frequency, respectively)
    % - FF       : [1 x N] time variant FF
     
    % Outputs:
    % - u      : control variable
    % - y      : criterion output
    % - dy     : criterion gradient estimate
    % - ddy    : hessian estimate
    % - ddy_inv: hessian inverse estimate

    if (length(A) ~= length(f)) || (length(f) ~= size(ddy0, 1))
        fprintf('Dimensions do not match\n.')
        return       
    end    
  
    % Data structures
    n              = length(f);
    u_hat          = zeros(n, N+1);
    u              = [u0, zeros(n, N)];
    y              = zeros(1, N);
    y_hat          = zeros(1, N);
    dy             = zeros(n, N);
    hpf            = zeros(1, N);
    zeta           = zeros(n, N);
    lpf_grad       = zeros(n, N);
    sigma          = zeros(n, n, N);
    lpf_hessian    = zeros(n, n, N);
    ddy            = zeros(n, n, N);
    ddy_inv        = zeros(n, n, N+1); % output of ricatti filter
    ddy_inv(:,:,1) = ddy0;
    
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

    % N(t) - Hessian estimate
    Nhess = cell(n,n);
    for i = 1:n
        for j = 1:n
            if i == j
                Nhess{i,j} = @(t) 16 / (rad2deg(A(i))^2) * (sin(2*pi*f(i)*t)^2 - 0.5);
            else
                Nhess{i,j} = @(t) 4 / (rad2deg(A(i)) * rad2deg(A(j))) * sin(2*pi*f(i)*t) * sin(2*pi*f(j)*t);
            end
        end
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

        % HPF
        if i >= bworder+1
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
            % Gradient
            zeta(:, i) = hpf(i) .* sin(2*pi*f*t);
            
            % LPF - Gradient
            if i >= M
                lpf_grad(:, i) = 1/al(end) * (-al(1:end-1) * lpf_grad(:, i-M+1:i-1)' + bl * zeta(:, i-M+1:i)');           
            else
                for j = 1:i
                    lpf_grad(:, i) = lpf_grad(:, i) + bl(end-j+1) .* zeta(:, i-j+1);
                end
        
                for j = 2:i
                    lpf_grad(:, i) = lpf_grad(:, i) - al(end-j+1) .* lpf_grad(:, i-j+1);
                end
                
                lpf_grad(:, i) = 1/al(end) .* lpf_grad(:, i);
            end
            
            dy(:, i) = lpf_grad(:, i);
            
            % Hessian 
            for j = 1:n
                for k = 1:n
                    sigma(j,k,i) = hpf(i) * Nhess{j,k}(t); % \hat{H} = N(t)y_hp
                end
            end

            % LPF - Hessian
            if i >= M
                for j = 1:M
                    lpf_hessian(:, :, i) = lpf_hessian(:, :, i) + bl(end-j+1).*sigma(:, :, i-j+1);
                end
        
                for j = 2:M
                    lpf_hessian(:, :, i) = lpf_hessian(:, :, i) - al(end-j+1).*lpf_hessian(:, :, i-j+1);
                end
                
                lpf_hessian(:, :, i) = 1/al(end) * lpf_hessian(:, :, i);
            else
                for j = 1:i
                    lpf_hessian(:, :, i) = lpf_hessian(:, :, i) + bl(end-j+1).*sigma(:, :, i-j+1);
                end
        
                for j = 2:i
                    lpf_hessian(:, :, i) = lpf_hessian(:, :, i) - al(end-j+1)*lpf_hessian(:, :, i-j+1);
                end
                
                lpf_hessian(:, :, i) = 1/al(end) .* lpf_hessian(:, :, i);
            end
    
            ddy(:, :, i) = lpf_hessian(:, :, i);
        
        else
            % Gradient
            dy(:, i) = hpf(i) * sin(2*pi*f*t);
            
            % Hessian
            for j = 1:n
                for k = 1:n
                    ddy(j,k,i) = hpf(i) * Nhess{j,k}(t); % \hat{H} = N(t)y_hp
                end
            end            
        end       
        
        % Euler discretization of Ricatti equation    
        ddy_inv(:, :, i+1) = ddy_inv(:, :, i) + dt * (wric * ddy_inv(:, :, i) - ...
            wric * ddy_inv(:, :, i) * ddy(:, :, i) * ddy_inv(:, :, i)); 

        % Parameter estimate - Single Integrator
        u_hat(:, i+1) = u_hat(:, i) - dt * K * ddy_inv(:, :, i+1) * dy(:, i); 
        
        % Add dither
        u(:, i+1) = FF(i) + u_hat(:, i+1) + A .* sin(2*pi*f*t);
        
        % Error condition
        if any(u(:, i+1) > pi) || any(u(:, i+1) < -pi) 
            break
        end
    end
    toc
end