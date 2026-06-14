%% MAIN

clc;
clear;
close all;
echo off; % Used to avoid unwanted warnings and other console stuff

addpath("functions\")
addpath("functions\lineoptifuncs\")

% Run Params
trackplot=false;

par = carParams();
n_var = 200;          % Number of Design Variables for Interpolation
car_margin = 0.5;    % (0.5m margin)

% Generate Track
track = genTrack(trackplot);

% Set up for optimizer
lineopti.s_full = [0; cumsum(track.vecmag(1:end-1))];        % Cumulative distance
lineopti.s_ctrl = linspace(0, lineopti.s_full(end), n_var)';

lineopti.w_left_ctrl = interp1(lineopti.s_full, track.w(:,1), lineopti.s_ctrl);
lineopti.w_right_ctrl = interp1(lineopti.s_full, track.w(:,2), lineopti.s_ctrl);

% fmincon lower and upper bounds (alpha is positive to the left)

lb = -lineopti.w_right_ctrl + car_margin; 
ub =  lineopti.w_left_ctrl  - car_margin;

% Initial guess (start exactly on the centerline, so alpha = 0)
lineopti.alpha_guess = zeros(n_var, 1);

% figure;
% plot(lineopti.s_full,track.w(:,1))
% hold on
% plot(lineopti.s_ctrl,lineopti.w_left_ctrl)

%%% OPTIMIZATION ROUTINE

% Define fmincon options
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 2000000, ... 
    'MaxIterations', 20, ...           
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

weight_length = 0.001;
objectiveFcn = @(alpha) calcCurvatureCost(alpha, lineopti.s_ctrl, lineopti.s_full, track, weight_length);

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
lineopti.alpha_opti_full = makima(lineopti.s_ctrl, lineopti.alpha_opt, lineopti.s_full);
bdeg = 3;
bknots = augknt(lineopti.s_ctrl, bdeg+1);
b_spline_curve = spmak(bknots, lineopti.alpha_opt');
lineopti.alpha_opti_full = fnval(b_spline_curve, lineopti.s_full);
lineopti.optimized = track.m + track.vecleft ./ track.vecmag .* lineopti.alpha_opti_full;

[laptime, vprof] = calcTimeAndVelocity(lineopti.alpha_opt, lineopti.s_ctrl, lineopti.s_full, track, par, 'bspline');

%%% Calculate Acceleration Profile
% Acceleration = v * dv/ds
dvds = gradient(vprof, lineopti.s_full);
accel = vprof .* dvds; 

%%% Calculate Acceleration Profile
% Acceleration = v * dv/ds
dvds = gradient(vprof, lineopti.s_full);
accel = vprof .* dvds; 

%%% Define Color Limits Explicitly
minvel = min(vprof);
maxvel = max(vprof);
minacc = min(accel);
maxacc = max(accel);

%% --- Plot 1: Optimized Line (Velocity) ---
figure('Name', 'Optimized Racing Line - Velocity', 'Color', 'w', 'Position', [100, 100, 900, 700]);
hold on; grid on; axis equal;

% 1. Fill the track surface
trackX = [track.l(:,1); flipud(track.r(:,1))];
trackY = [track.l(:,2); flipud(track.r(:,2))];
fill(trackX, trackY, [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');

% 2. Plot Track Boundaries & Centerline
plot(track.l(:,1), track.l(:,2), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Track Limits');
plot(track.r(:,1), track.r(:,2), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(track.m(:,1), track.m(:,2), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'Centerline');

% 3. Plot Optimized Line (Color-coded by Velocity)
x = lineopti.optimized(:,1)';
y = lineopti.optimized(:,2)';
z = zeros(size(x));
col_vel = vprof(:)';

surface([x;x], [y;y], [z;z], [col_vel;col_vel], ...
    'FaceColor', 'no', ...
    'EdgeColor', 'interp', ...
    'LineWidth', 2.5);

plot(NaN, NaN, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Optimized Line');

% 4. Formatting
colormap(gca, parula); 
clim([minvel, maxvel]); % Apply strict min/max scaling for velocity
c = colorbar;
c.Label.String = 'Velocity [m/s]';
c.Label.FontSize = 11;
c.Label.FontWeight = 'bold';

xlabel('X Position [m]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y Position [m]', 'FontSize', 11, 'FontWeight', 'bold');
title('Optimal Racing Line - Velocity Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);


%% --- Plot 2: Optimized Line (Acceleration / Deceleration) ---
figure('Name', 'Optimized Racing Line - Acceleration', 'Color', 'w', 'Position', [150, 150, 900, 700]);
hold on; grid on; axis equal;

% 1. Fill the track surface
fill(trackX, trackY, [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');

% 2. Plot Track Boundaries & Centerline
plot(track.l(:,1), track.l(:,2), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Track Limits');
plot(track.r(:,1), track.r(:,2), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(track.m(:,1), track.m(:,2), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'Centerline');

% 3. Plot Optimized Line (Color-coded by Acceleration)
col_acc = accel(:)';

surface([x;x], [y;y], [z;z], [col_acc;col_acc], ...
    'FaceColor', 'no', ...
    'EdgeColor', 'interp', ...
    'LineWidth', 2.5);

plot(NaN, NaN, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Optimized Line');

% 4. Formatting
colormap(gca, turbo); 

% Apply strict min/max scaling for acceleration 
clim([minacc, maxacc]); 
% NOTE: If you notice the 0-point (coasting) isn't dead-center in the color scale 
% because braking is stronger than acceleration, you can change this line to: 
% clim([-max(abs(accel)), max(abs(accel))]); to force symmetry.

c2 = colorbar;
c2.Label.String = 'Acceleration [m/s^2]';
c2.Label.FontSize = 11;
c2.Label.FontWeight = 'bold';

xlabel('X Position [m]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y Position [m]', 'FontSize', 11, 'FontWeight', 'bold');
title('Optimal Racing Line - Acceleration & Braking Zones', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);


%% --- Plot 3: 2D Velocity & Acceleration Graphs ---
figure('Name', 'Telemetry Profiles', 'Color', 'w', 'Position', [200, 200, 800, 500]);

% Velocity Subplot
subplot(2,1,1);
hold on; grid on; box on;
area(lineopti.s_full, vprof, 'FaceColor', [0.85 0.9 0.95], 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(lineopti.s_full, vprof, '-', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
xlim([0, lineopti.s_full(end)]);
ylim([max(0, minvel - 2), maxvel * 1.1]); % Slightly adjusted to use minvel
ylabel('Velocity [m/s]', 'FontSize', 11, 'FontWeight', 'bold');
title('Telemetry along Racing Line', 'FontSize', 14, 'FontWeight', 'bold');

% Acceleration Subplot
subplot(2,1,2);
hold on; grid on; box on;
plot(lineopti.s_full, accel, '-', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
yline(0, 'k--', 'LineWidth', 1); % 0 line to distinguish accel/brake
xlim([0, lineopti.s_full(end)]);
ylim([minacc * 1.1, maxacc * 1.1]); % Set Y-limits based on actual min/max
xlabel('Track Distance [m]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Accel [m/s^2]', 'FontSize', 11, 'FontWeight', 'bold');
%%% Calculate Line Parameters

lineopti.dx = gradient(lineopti.optimized(:,1));
lineopti.dy = gradient(lineopti.optimized(:,2));
lineopti.ddx = gradient(lineopti.dx);
lineopti.ddy = gradient(lineopti.dy);
lineopti.kappasquare = sum(((lineopti.dx .* lineopti.ddy - lineopti.dy .* lineopti.ddx) ./ ((lineopti.dx.^2 + lineopti.dy.^2).^(3/2))).^2);
lineopti.length = sum(sqrt(lineopti.dx.^2 + lineopti.dy.^2));

fprintf("Track Midline Square Curvature: %e\n", track.kappasquare)
fprintf("Track Midline Length [m]: %e\n", track.length)
fprintf("Optimized Line Square Curvature: %e\n", lineopti.kappasquare)
fprintf("Line Length [m]: %e\n", lineopti.length)

[laptime,vprof] = calcTimeAndVelocity(lineopti.alpha_opt, lineopti.s_ctrl, lineopti.s_full, track, par, 'bspline');
fprintf("Lap Time [s]: %e\n", laptime)

figure;
plot(lineopti.s_full,vprof)