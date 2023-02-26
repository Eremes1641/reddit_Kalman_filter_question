function [stateDt,simOut] = plant_pendulum(t,state,parm)
    %% this plant is a pendulum without coulomb fricion and drived by voltage
    %% declare
    % plant
    a = parm.a;  % viscous friction [N-m / rad/s]
    b = parm.b;    % gain [N-m/volt]
    g = parm.g;
    l = parm.l; % length [m]

    % u
    t0 = 0.1;
    t1 = 0.5;

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
    d_gravity = g/l*sin(x1);

    %% plant
    x1Dt = x2;
    x2Dt = -a*x2 + b*u - d_gravity;

    %% return
    stateDt = [x1Dt;
        x2Dt];
    simOut.t = t;
    simOut.u = u;       % input [volt]
    simOut.x1 = x1;     % pos [rad]
    simOut.x2 = x2;     % vel [rad/s]
    simOut.acc = x2Dt;  % acceleration [rad/s/s]
    simOut.d = d_gravity; % disturbance [volt]
end