%% MAIN
clc;
clear;
close all;
echo off; 
addpath("functions\")
addpath("functions\lineoptifuncs\")

% Run params
trackplot = false;
par = carParams();
splineType = 'makima';

%setup parameters individually
n_var_pso  = 50;   % for the pso we used a coarser grid, the get a bit more general solution
n_var_fmin = 350;  
car_margin = 0.5;  

% Generate Track
track = genTrack(trackplot);

% baseline distance
lineopti.s_full = [0; cumsum(track.vecmag(1:end-1))]; 

corner_weight = 40;    
smoothing_window = 30; 

% pso grid
lineopti.s_ctrl_pso = genAdaptiveNodes(lineopti.s_full, track.m, n_var_pso, corner_weight, smoothing_window);
lineopti.w_left_pso  = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl_pso, 'makima');
lineopti.w_right_pso = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl_pso, 'makima');
lb_pso = -lineopti.w_right_pso + car_margin; 
ub_pso =  lineopti.w_left_pso  - car_margin;
Aeq_pso = zeros(2, n_var_pso);
Aeq_pso(1, 1) = 1; Aeq_pso(1, end) = -1;
Aeq_pso(2, 1) = -1; Aeq_pso(2, 2) = 1; Aeq_pso(2, end-1) = 1; Aeq_pso(2, end) = -1;

% fmincon grid
lineopti.s_ctrl_fmin = genAdaptiveNodes(lineopti.s_full, track.m, n_var_fmin, corner_weight, smoothing_window);
lineopti.w_left_fmin  = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl_fmin, 'makima');
lineopti.w_right_fmin = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl_fmin, 'makima');

lb_fmin = -lineopti.w_right_fmin + car_margin; 
ub_fmin =  lineopti.w_left_fmin  - car_margin;

Aeq_fmin = zeros(2, n_var_fmin);
Aeq_fmin(1, 1) = 1; Aeq_fmin(1, end) = -1;
Aeq_fmin(2, 1) = -1; Aeq_fmin(2, 2) = 1; Aeq_fmin(2, end-1) = 1; Aeq_fmin(2, end) = -1;

beq = zeros(2, 1); 

% perform PSO
fprintf('\n--- Starting Phase 1: PSO (Global Search on %d nodes) ---\n', n_var_pso);
n_particles = 100;     
n_iterations = 600;   
w = 0.7; cp = 1.5; cg = 1.5; 
verbose = true;       

penalty_weight = 1e5;

objectiveFcnPSO = @(alpha) calcLapTimeCostDetailFriction(alpha, lineopti.s_ctrl_pso, lineopti.s_full, track, par, splineType) ...
                           + penalty_weight * sum((Aeq_pso * alpha - beq).^2);

[alpha_pso_rough, best_cost_pso, gs_hist, eval_hist] = pso(...
    objectiveFcnPSO, n_var_pso, lb_pso, ub_pso, n_particles, n_iterations, w, cp, cg, verbose);

fprintf('\nPhase 1 Complete. Bridging grids...\n');
alpha_fmin_guess = makima(lineopti.s_ctrl_pso, alpha_pso_rough, lineopti.s_ctrl_fmin);
alpha_fmin_guess = max(min(alpha_fmin_guess, ub_fmin), lb_fmin); % Safety clamp
%fmincon phase
fprintf('\n--- Starting Phase 2: fmincon (Local Smoothing on %d nodes) ---\n', n_var_fmin);

objectiveFcnExact = @(alpha) calcLapTimeCostDetailFriction(alpha, lineopti.s_ctrl_fmin, lineopti.s_full, track, par, splineType);

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...       
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 20000000, ...
    'MaxIterations', 2000, ...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

lineopti.alpha_opt = fmincon(objectiveFcnExact, alpha_fmin_guess, [], [], Aeq_fmin, beq, lb_fmin, ub_fmin, [], options);

fprintf('\n--- Hybrid Optimization Complete ---\n');

%post processing (plotting not up to date)
if isequal(splineType,'makima')
    lineopti.alpha_opti_full = makima(lineopti.s_ctrl_fmin, lineopti.alpha_opt, lineopti.s_full);
elseif isequal(splineType,'bspline')
    bdeg = 3;
    bknots = augknt(lineopti.s_ctrl_fmin, bdeg+1);
    b_spline_curve = spmak(bknots, lineopti.alpha_opt');
    lineopti.alpha_opti_full = fnval(b_spline_curve, lineopti.s_full);
end

lineopti.optimized = track.m + track.vecleft ./ track.vecmag .* lineopti.alpha_opti_full;

figure;
plot(track.m(:,1), track.m(:,2), lineopti.optimized(:,1), lineopti.optimized(:,2), track.l(:,1), track.l(:,2), track.r(:,1), track.r(:,2))
legend("Centerline", "Optimized Line", "Left Bound", "Right Bound")

%%% Calculate Line Parameters
lineopti.dx = gradient(lineopti.optimized(:,1));
lineopti.dy = gradient(lineopti.optimized(:,2));
lineopti.ddx = gradient(lineopti.dx);
lineopti.ddy = gradient(lineopti.dy);
lineopti.kappasquare = sum(((lineopti.dx .* lineopti.ddy - lineopti.dy .* lineopti.ddx) ./ ((lineopti.dx.^2 + lineopti.dy.^2).^(3/2))).^2);
lineopti.ds = sqrt(lineopti.dx.^2 + lineopti.dy.^2);
lineopti.length = sum(lineopti.ds);

fprintf("Track Midline Square Curvature: %e\n", track.kappasquare)
fprintf("Track Midline Length [m]: %e\n", track.length)
fprintf("Optimized Line Square Curvature: %e\n", lineopti.kappasquare)
fprintf("Line Length [m]: %e\n", lineopti.length)
[laptime, vprof] = calcTimeAndVelocity(lineopti.alpha_opt, lineopti.s_ctrl_fmin, lineopti.s_full, track, par, splineType);

fprintf("Lap Time [s]: %e\n", laptime)

figure;
plot(lineopti.s_full, vprof)