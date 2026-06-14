%%% SET CAR PARAMETERS

function par = carParams()

    % Simple Params
    % par.Vmax = 320/3.6; % m/s
    par.m = 798; % kg
    par.P_engine = 735500; % W, 1000hp
    % par.F_engine = 5500; % N
    par.F_engine_max = 10000;
    par.R_min = 5.0;
    
    % Second Order Params - Required for more advanced model
    par.CLA = 4.8;
    par.CDA = 1.2;
    par.rho = 1.225;
    % par.accelmax = 1.5; % g
    % par.decelmax = 5; % g
    % par.latmax = 5; % g

    par.Vmax = (par.P_engine/(par.CDA*0.5*par.rho))^(1/3)

end