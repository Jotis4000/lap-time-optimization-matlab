function cost_time = calcLapTimePSO(alpha_ctrl, s_ctrl, s_full, track, par, splineType)

    % size(alpha_ctrl)
    % size(s_ctrl)
    % size(s_full)
    % alpha_ctrl = alpha_ctrl';

    % Interpolate to the full 900 points
    if isequal(splineType,'makima')
        alpha_full = makima(s_ctrl, alpha_ctrl, s_full);
    elseif isequal(splineType,'bspline')
        alpha_ctrl = alpha_ctrl';
        bdeg = 3;
        bknots = augknt(s_ctrl,bdeg+1);
        b_spline_curve = spmak(bknots, alpha_ctrl');
        alpha_full = fnval(b_spline_curve, s_full);
    end
    
    % track limits
    car_margin = 0.5;
    c_left  = alpha_full - track.w(:,1) + car_margin;
    c_right = -alpha_full - track.w(:,2) + car_margin;
    
    % Find the maximum distance the car went out of bounds anywhere on track
    max_violation = max([c_left; c_right; 0]); 
    
    % Calculate Geometry
    nx = track.vecleft(:,1) ./ track.vecmag;
    ny = track.vecleft(:,2) ./ track.vecmag;
    X_race = track.m(:,1) + alpha_full .* nx;
    Y_race = track.m(:,2) + alpha_full .* ny;

    dx = gradient(X_race);
    dy = gradient(Y_race);
    ddx = gradient(dx);
    ddy = gradient(dy);
    
    kappa = (dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2));
    % ds_race = sqrt(dx.^2 + dy.^2); 
    
    %steering lim
    kappa_max = 1 / par.R_min;
    steering_violation = max(max(abs(kappa)) - kappa_max, 0);
    % Run physics model
    actual_lap_time = calcLapTimeCostDetail(alpha_ctrl, s_ctrl, s_full, track, par, splineType);
    
    % Penalties
    if max_violation > 0
        cost_time = actual_lap_time + (20.0 * max_violation);
    else if alpha_full(1) ~= alpha_full(end)
        cost_time = actual_lap_time + (20.0 * max_violation);
    else
        cost_time = actual_lap_time;
    end
    
    if steering_violation > 0
        cost_time = cost_time + (50.0 * steering_violation);
    end
end