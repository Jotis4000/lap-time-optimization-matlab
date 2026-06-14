function track = genTrack(plt)

    if nargin < 1
        plt = true;
    end
    %% Parameters
    addpath("tracks\")

    FileName = "Budapest.csv";
    datat = xlsread(FileName); % x_m, y_m, w_right, w_left

    track.m = [datat(:,1) datat(:,2)]; % X, Y
    track.w = [datat(:,4) datat(:,3)]; % Left, Right

    % generate track bounds: iterate through data set
    track.vec = [circshift(track.m(:,1), [-1 0])-track.m(:,1) circshift(track.m(:,2), [-1 0])-track.m(:,2)];
    track.vecmag = vecnorm(track.vec,2,2);
    track.vecleft = [-track.vec(:,2), track.vec(:,1)];
    track.vecright = [track.vec(:,2), -track.vec(:,1)];
    track.l = track.m+track.vecleft./track.vecmag.*track.w(:,1);
    track.r = track.m+track.vecright./track.vecmag.*track.w(:,2);
    track.length = sum(track.vecmag);

    if plt
        plot(track.m(:,1),track.m(:,2),track.l(:,1),track.l(:,2),track.r(:,1),track.r(:,2))
        legend("m","l","r")
    end

    % Calc midline curvature
    dx = gradient(track.m(:,1));
    dy = gradient(track.m(:,2));
    ddx = gradient(dx);
    ddy = gradient(dy);
    track.kappasquare = sum(((dx .* ddy - dy .* ddx) ./ ((dx.^2 + dy.^2).^(3/2))).^2);

    track.mu = 1.7; % Ground friction coefficient

end

