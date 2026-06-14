%% MAIN
clc;
clear;
close all;
echo off; % Used to avoid unwanted warnings and other console stuff

addpath("functions\")

% Run Params
trackplot=false;

par = carParams();
n_var = 300;          % Number of Design Variables for Interpolation
car_margin = 0.5;    % Car half-width margin (0.5m margin)

guess = 'geom'; % 'center', 'geom'
splineType = 'makima'; % 'makima', 'bspline'

% Generate Track
track = genTrack(trackplot);

% Set up for optimizer
lineopti.s_full = [0; cumsum(track.vecmag(1:end-1))];        % Cumulative distance
lineopti.s_ctrl = linspace(0, lineopti.s_full(end), n_var)';

corner_weight = 20;    % How strongly corners pull nodes (Multiplier)
smoothing_window = 50; % How wide the apex cluster is
lineopti.s_ctrl = genAdaptiveNodes(lineopti.s_full, track.m, n_var, corner_weight, smoothing_window);

lineopti.w_left_ctrl = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl);
lineopti.w_right_ctrl = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl);

% fmincon lower and upper bounds (alpha is positive to the left)

if isequal(splineType,'makima')
    bmargin = 0;
elseif isequal(splineType,'bspline')
    bmargin = 10;
end
lb = -lineopti.w_right_ctrl + car_margin-bmargin; 
ub =  lineopti.w_left_ctrl  - car_margin+bmargin;

%%% OPTIMIZATION ROUTINE

% Initial guess (start exactly on the centerline, so alpha = 0)
lineopti.alpha_guess = zeros(n_var, 1);

if isequal(guess,'center')
    lineopti.alpha_guess = zeros(n_var, 1);
elseif isequal(guess,'geom')
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 2000000, ... 
        'MaxIterations', 1000, ...           
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6);

    weight_length = 0;
    objectiveFcn = @(alpha) calcCurvatureCost(alpha, lineopti.s_ctrl, lineopti.s_full, track, weight_length);

    Aeq = zeros(2, n_var);
    beq = zeros(2, 1); % The right side of the equations (both equal 0)
    
    Aeq(1, 1)   = 1;
    Aeq(1, end) = -1;
    
    Aeq(2, 1)     = -1;
    Aeq(2, 2)     =  1;
    Aeq(2, end-1) =  1;
    Aeq(2, end)   = -1;
    
    % x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
    lineopti.alpha_opt = fmincon(objectiveFcn, lineopti.alpha_guess, [], [], Aeq, beq, lb, ub, [], options);
    lineopti.alpha_guess = lineopti.alpha_opt;
end

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 20000000, ... 
    'MaxIterations', 2000, ...           
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6); %, 'FiniteDifferenceType', 'central');

objectiveFcn = @(alpha) calcLapTimeCostDetailFriction(alpha, lineopti.s_ctrl, lineopti.s_full, track, par, splineType);

% Equality constraints
% First and last point need to have the same position and direction
Aeq = zeros(2, n_var);
beq = zeros(2, 1); % The right side of the equations (both equal 0)

Aeq(1, 1)   = 1;
Aeq(1, end) = -1;

Aeq(2, 1)     = -1;
Aeq(2, 2)     =  1;
Aeq(2, end-1) =  1;
Aeq(2, end)   = -1;

% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
lineopti.alpha_opt = fmincon(objectiveFcn, lineopti.alpha_guess, [], [], Aeq, beq, lb, ub, [], options);

%%% Generate optimized line and plot result

if isequal(splineType,'makima')
    lineopti.alpha_opti_full = makima(lineopti.s_ctrl, lineopti.alpha_opt, lineopti.s_full);
elseif isequal(splineType,'bspline')
    bdeg = 3;
    bknots = augknt(lineopti.s_ctrl,bdeg+1);
    b_spline_curve = spmak(bknots, lineopti.alpha_opt');
    lineopti.alpha_opti_full = fnval(b_spline_curve, lineopti.s_full);
end

lineopti.optimized = track.m+track.vecleft./track.vecmag.*lineopti.alpha_opti_full;

figure;
plot(track.m(:,1),track.m(:,2),lineopti.optimized(:,1),lineopti.optimized(:,2),track.l(:,1),track.l(:,2),track.r(:,1),track.r(:,2))
legend("Centerline","Optimized Line","Left Bound","Right Bound")

%%% Calculate Line Parameters

lineopti.dx = gradient(lineopti.optimized(:,1));
lineopti.dy = gradient(lineopti.optimized(:,2));
lineopti.ddx = gradient(lineopti.dx);
lineopti.ddy = gradient(lineopti.dy);
lineopti.kappasquare = sum(((lineopti.dx .* lineopti.ddy - lineopti.dy .* lineopti.ddx) ./ ((lineopti.dx.^2 + lineopti.dy.^2).^(3/2))).^2);
lineopti.ds = sqrt(lineopti.dx.^2 + lineopti.dy.^2);
lineopti.length = sum(lineopti.ds);
% lineopti.laptime = sum(ds./v_profile);

fprintf("Track Midline Square Curvature: %e\n", track.kappasquare)
fprintf("Track Midline Length [m]: %e\n", track.length)
fprintf("Optimized Line Square Curvature: %e\n", lineopti.kappasquare)
fprintf("Line Length [m]: %e\n", lineopti.length)
% lpha_ctrl, s_ctrl, s_full, track, par
[laptime, vprof] = calcTimeAndVelocity(lineopti.alpha_opt, lineopti.s_ctrl, lineopti.s_full, track, par, splineType);
fprintf("Lap Time [s]: %e\n", laptime)

figure;
plot(lineopti.s_full,vprof)