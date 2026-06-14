function cost = calcCurvatureCost(alpha_ctrl, s_ctrl, s_full, track, weight_length)
    alpha_full = makima(s_ctrl, alpha_ctrl, s_full);
    % bdeg = 3;
    % bknots = augknt(s_ctrl,bdeg+1);
    % b_spline_curve = spmak(bknots, alpha_ctrl');
    % alpha_full = fnval(b_spline_curve, s_full);
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
    ds_race = sqrt(dx.^2 + dy.^2); 
    
    % Caclulate costs
    cost_curvature = sum((kappa.^2) .* ds_race); 
    cost_length    = sum(ds_race);               % Total distance traveled
    
    cost = (cost_curvature + (weight_length * cost_length)); % 1000 scaling to reduce function operations
    
end