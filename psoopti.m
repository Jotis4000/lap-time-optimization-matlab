%% MAIN
clc;
clear;
close all;
echo off; 

addpath("functions\")
addpath("functions\lineoptifuncs\")

% run params
trackplot = false;
par = carParams();
n_var = 160;          
car_margin = 0.5;  %carhalf margin

%Generate track
track = genTrack(trackplot);

lineopti.s_full = [0; cumsum(track.vecmag(1:end-1))];       
lineopti.s_ctrl = linspace(0, lineopti.s_full(end), n_var)';
lineopti.w_left_ctrl = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl);
lineopti.w_right_ctrl = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl);

% Lower and upper bounds
lb = -lineopti.w_right_ctrl + car_margin; 
ub =  lineopti.w_left_ctrl  - car_margin;

% eq constraits
Aeq = zeros(2, n_var);
beq = zeros(2, 1); 
Aeq(1, 1)   = 1;
Aeq(1, end) = -1;
Aeq(2, 1)     = -1;
Aeq(2, 2)     =  1;
Aeq(2, end-1) =  1;
Aeq(2, end)   = -1;

%%% PSO opt

%hyperparameters
n_particles = 80;     % swarmsize
n_iterations = 200;   % maxiter
w = 0.7;              % inertia
cp = 1.5;             % cognitive param
cg = 1.5;             % social param
verbose = true; 
% objective func with constraints and penalty
penalty_weight = 1e6;
objectiveFcnPSO = @(alpha) calcLapTimeCostDetail(alpha, lineopti.s_ctrl, lineopti.s_full, track, par, 'makima') + penalty_weight * sum((Aeq * alpha - beq).^2);
% run pso
[lineopti.alpha_opt, best_cost, gs_hist, eval_hist] = pso(objectiveFcnPSO, n_var, lb, ub, n_particles, n_iterations, w, cp, cg, verbose);

% Generate optimized line and plot result
lineopti.alpha_opti_full = makima(lineopti.s_ctrl, lineopti.alpha_opt, lineopti.s_full);
lineopti.optimized = track.m + track.vecleft ./ track.vecmag .* lineopti.alpha_opti_full;

figure('Name', 'Track Optimization Result');
plot(track.m(:,1), track.m(:,2), '--k', ...
     lineopti.optimized(:,1), lineopti.optimized(:,2), '-b', ...
     track.l(:,1), track.l(:,2), '-r', ...
     track.r(:,1), track.r(:,2), '-r', 'LineWidth', 1.5)
legend("Centerline", "Optimized Line", "Track Bounds", "", "Location", "best")
axis equal;
title('Optimized Racing Line (PSO)');

%%% Calculate Line Parameters
lineopti.dx = gradient(lineopti.optimized(:,1));
lineopti.dy = gradient(lineopti.optimized(:,2));
lineopti.ddx = gradient(lineopti.dx);
lineopti.ddy = gradient(lineopti.dy);
lineopti.kappasquare = sum(((lineopti.dx .* lineopti.ddy - lineopti.dy .* lineopti.ddx) ./ ((lineopti.dx.^2 + lineopti.dy.^2).^(3/2))).^2);
lineopti.ds = sqrt(lineopti.dx.^2 + lineopti.dy.^2);
lineopti.length = sum(lineopti.ds);

fprintf("\n--- FINAL METRICS ---\n")
fprintf("Track Midline Square Curvature: %e\n", track.kappasquare)
fprintf("Track Midline Length [m]: %e\n", track.length)
fprintf("Optimized Line Square Curvature: %e\n", lineopti.kappasquare)
fprintf("Line Length [m]: %e\n", lineopti.length)

% Calculate final lap time
[laptime, vprof] = calcLapTimeCostDetail(lineopti.alpha_opt, lineopti.s_ctrl, lineopti.s_full, track, par);
actual_laptime = calcLapTimeCostDetail(lineopti.alpha_opt, lineopti.s_ctrl, lineopti.s_full, track, par);
fprintf("Estimated Lap Time Cost: %e\n", actual_laptime)

% Plot PSO Convergence
figure('Name', 'PSO Convergence');
plot(1:n_iterations, eval_hist, 'LineWidth', 2);
title('PSO Convergence History');
xlabel('Iteration');
ylabel('Penalized Cost');
grid on;