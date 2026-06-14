function s_ctrl = genAdaptiveNodes(s_full, track_m, n_var, corner_weight, smoothing_window)
    
    % Calculate raw curvature 
    dx = gradient(track_m(:,1));
    dy = gradient(track_m(:,2));
    ddx = gradient(dx);
    ddy = gradient(dy);
    kappa_raw = (dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2));

    % smooth
    kappa_smooth = smoothdata(abs(kappa_raw), 'gaussian', smoothing_window); 

    % normalize
    kappa_norm = kappa_smooth / max(kappa_smooth);
    kappa_focused = kappa_norm .^ 2; 

    % node density
    node_density = 1 + (corner_weight .* kappa_focused);
    cum_density = cumtrapz(s_full, node_density);

    % sample
    cum_density_samples = linspace(0, cum_density(end), n_var)';
    
    % map into domain
    s_ctrl = interp1(cum_density, s_full, cum_density_samples, 'pchip');
end