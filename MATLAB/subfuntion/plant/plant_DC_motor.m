function [stateDt,simOut] = plant_DC_motor(t,state,parm)
    %% this plant is a DC motor with coulomb fricion drived by voltage
    %% declare
    % plant
    a = parm.a;  % viscous friction [N-m / rad/s]
    b = parm.b;  % gain [N-m/volt]
    d_Coulomb_coeff = parm.d_Coulomb_coeff;         % Coulomb friction coeff [volt]
    d_Coulomb_threshold = parm.d_Coulomb_threshold; % Coulomb friction thrshould [rad/s]

    % u
    t0 = 0.1;
    t1 = 1;

    %% substitute
    % plant state
    x1 = state(1);  % pos [rad]
    x2 = state(2);  % vel [rad/s]

    %% u (volt)
    if t >= t0 && t < t1
        u = 10;
    else
        u = 0;
    end

    %% disturbance, Coulomb friction [volt]
    d_Coulomb = d_Coulomb_coeff*tanh(x2/d_Coulomb_threshold);

    %% plant
    x1Dt = x2;
    x2Dt = -a*x2 + b*u - b*d_Coulomb;

    %% return
    stateDt = [x1Dt;
        x2Dt];
    simOut.t = t;
    simOut.u = u;       % input [volt]
    simOut.x1 = x1;     % pos [rad]
    simOut.x2 = x2;     % vel [rad/s]
    simOut.acc = x2Dt;  % acceleration [rad/s/s]
    simOut.d = d_Coulomb; % disturbance [volt]
end