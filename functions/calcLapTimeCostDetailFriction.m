function cost = calcLapTimeCostDetailFriction(alpha_ctrl, s_ctrl, s_full, track, par, splineType)
    
    % addpath("functions\") this mfer <----
    N = length(s_full);
    if isequal(splineType,'makima')
        alpha_full = makima(s_ctrl, alpha_ctrl, s_full);
    elseif isequal(splineType,'bspline')
        bdeg = 3;
        bknots = augknt(s_ctrl,bdeg+1);
        b_spline_curve = spmak(bknots, alpha_ctrl');
        alpha_full = fnval(b_spline_curve, s_full);
    end
    
    % Calculate line parameters the same as in the geometriic optimization
    nx = track.vecleft(:,1) ./ track.vecmag;
    ny = track.vecleft(:,2) ./ track.vecmag;
    X_race = track.m(:,1) + alpha_full .* nx;
    Y_race = track.m(:,2) + alpha_full .* ny;
    
    dx = gradient(X_race);
    dy = gradient(Y_race);
    ddx = gradient(dx);
    ddy = gradient(dy);
    
    % Calculate curvature
    kappa = (dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2));
    % Calculate length
    ds = sqrt(dx.^2 + dy.^2);
    g = 9.81;
    
    Lconst = 0.5 * par.rho * par.CLA; % Downforce const
    Dconst = 0.5 * par.rho * par.CDA; % Drag const
    
    %% Aero Grip Ceiling
    denominator = (par.m .* abs(kappa)) - (track.mu * Lconst);
    v2_profile = zeros(N, 1);
    
    for i = 1:N
        if denominator(i) <= 0
            % Downforce outscales cornering demand
            v2_profile(i) = par.Vmax^2;
        else
            % Calculate limit, cap at car top speed
            v2_limit = (track.mu * par.m * g) / denominator(i);
            v2_profile(i) = min(v2_limit, par.Vmax^2);
        end
    end
    
    %% lap stiching loop
    % We run the backward and forward passes twice to wrap the start/finish line
    for pass_iter = 1:2
        
        %% Backward Integration (brakeing)
        dist_N_to_1 = sqrt((X_race(1)-X_race(N))^2 + (Y_race(1)-Y_race(N))^2);
        
        v2_next_bridge = v2_profile(1);
        F_normal_bridge = par.m * g + (Lconst * v2_next_bridge);
        F_grip_max_bridge = track.mu * F_normal_bridge;
        
        % friction circle
        F_lat_bridge = par.m * abs(kappa(1)) * v2_next_bridge;
        F_lon_avail_bridge = sqrt(max(0, F_grip_max_bridge^2 - F_lat_bridge^2));
        F_brake_total_bridge = F_lon_avail_bridge + (Dconst * v2_next_bridge);
        a_brake_bridge = F_brake_total_bridge / par.m;
        delta_v2_bridge = 2 * a_brake_bridge * dist_N_to_1;
        
        v2_profile(N) = min(v2_profile(N), v2_next_bridge + delta_v2_bridge);
        
        % backward pass
        for i = (N-1):-1:1
            v2_next = v2_profile(i+1);
            F_normal = par.m * g + (Lconst * v2_next);
            F_grip_max = track.mu * F_normal;
            
            % friction circle
            F_lat = par.m * abs(kappa(i+1)) * v2_next;
            F_lon_avail = sqrt(max(0, F_grip_max^2 - F_lat^2));
            F_brake_total = F_lon_avail + (Dconst * v2_next);
            a_brake = F_brake_total / par.m;
            delta_v2 = 2 * a_brake * ds(i);
            v2_profile(i) = min(v2_profile(i), v2_next + delta_v2);
        end
        
        %% Forward Integration (accel)
        v2_prev_bridge = v2_profile(N);
        v_prev_bridge = sqrt(v2_prev_bridge); 
        F_downforce_bridge = Lconst * v2_prev_bridge;
        F_drag_bridge      = Dconst * v2_prev_bridge;
        F_normal_bridge = par.m * 9.81 + F_downforce_bridge;
        F_grip_max_bridge = track.mu * F_normal_bridge;
        
        % friction circle
        F_lat_bridge = par.m * abs(kappa(N)) * v2_prev_bridge;
        F_lon_avail_bridge = sqrt(max(0, F_grip_max_bridge^2 - F_lat_bridge^2));
        
        if v_prev_bridge < 1.0 
            F_engine_bridge = par.F_engine_max;
        else
            F_engine_bridge = min(par.P_engine / v_prev_bridge, par.F_engine_max);
        end
        
        F_thrust_net_bridge = F_engine_bridge - F_drag_bridge;
        F_accel_bridge = min(F_lon_avail_bridge, F_thrust_net_bridge);
        a_accel_bridge = F_accel_bridge / par.m;
        delta_v2_bridge = 2 * a_accel_bridge * dist_N_to_1;
        
        v2_profile(1) = min(v2_profile(1), v2_prev_bridge + delta_v2_bridge);
        
        % forward pass
        for i = 2:N
            v2_prev = v2_profile(i-1);
            v_prev = sqrt(v2_prev); 
            F_downforce = Lconst * v2_prev;
            F_drag      = Dconst * v2_prev;
            F_normal = par.m * 9.81 + F_downforce;
            F_grip_max = track.mu * F_normal;
            
            % friction circle
            F_lat = par.m * abs(kappa(i-1)) * v2_prev;
            F_lon_avail = sqrt(max(0, F_grip_max^2 - F_lat^2));
            
            if v_prev < 1.0 
                F_engine = par.F_engine_max;
            else
                F_engine = min(par.P_engine / v_prev, par.F_engine_max);
            end
            
            F_thrust_net = F_engine - F_drag;
            F_accel = min(F_lon_avail, F_thrust_net);
            a_accel = F_accel / par.m;
            delta_v2 = 2 * a_accel * ds(i-1);
            
            v2_profile(i) = min(v2_profile(i), v2_prev + delta_v2);
        end
        
    end
    
    v_profile = sqrt(v2_profile);
    laptime_cost = 10*sum(ds./v_profile);
    alpha_wiggle = diff(alpha_ctrl, 2);
    smoothness_penalty = sum(alpha_wiggle.^2);
    
    cost = laptime_cost + (0.005 * smoothness_penalty);
end