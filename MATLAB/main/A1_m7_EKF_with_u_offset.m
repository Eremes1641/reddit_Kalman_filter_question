clc; clear; close all;

%% declare
% simulation
Fs = 500;   % sampling freq [Hz]
dt = 1/Fs;  % sampling period [sec]

% plant
plant_parm.a = 4;  % viscous friction [N-m / rad/s]
plant_parm.b = 35; % gain [N-m/volt]
plant_parm.d_Coulomb_coeff = 5;       % Coulomb friction coeff [volt]
plant_parm.d_Coulomb_threshold = 3;   % Coulomb friction thrshould [rad/s]

% model
model_parm.a = 3;  % viscous friction [N-m / rad/s]
model_parm.b = 30; % gain [N-m/volt]
model_parm.d_Coulomb_coeff = 6;       % Coulomb friction coeff [volt]
model_parm.d_Coulomb_threshold = 2;   % Coulomb friction thrshould [rad/s]

% EKF parameter
Q = diag([1E3 1E3]); % u weighting, general disturbance weighting
R = 1E3;             % pos weighting

%% ODE
% simulation
time = 0:dt:1.5;
IC = [0 0]';
equ = @(t,IC)plant_DC_motor(t,IC,plant_parm);
[t_sim,X_sim] = ode45(equ, time, IC);

% extract data
simData = extract_sim_data(equ, t_sim, X_sim);
t = simData.t;
pos_true = simData.x1;
vel_true = simData.x2;
u_input = simData.u;
lenT = length(t);
% d = simData.d;
d_true = get_d_true(t_sim,X_sim,u_input,plant_parm, model_parm);

%% measurement data
% discrete position measurement
resolution = 1/100; % 100 CPR
pos_measured = quantize_signal(pos_true,resolution);

%% kalman filter
G = [0 0; 1 0; 0 1];

% init
I = eye(3);
X_prior = [0 0 0]';
P_prior = zeros(3);
Xest = zeros(lenT,3);
% run EKF
for i = 1:lenT
    % substitute
    y = pos_measured(i);              % measured pos
    u = u_input(i);                   % input
    h = plant_h(X_prior,model_parm);  % estimated pos
    H = phpX(X_prior,model_parm);     % measurement matrix
    F = pfpX(X_prior,u,dt,model_parm); % process matrix
    
    % measurement update
    K = P_prior*H'/(H*P_prior*H'+R);
    X_post = X_prior + K*(y-h);
    P_post = (I-K*H)*P_prior;
    
    % time update
    X_prior = plant_f(X_post,u,dt,model_parm);
    P_prior = F*P_post*F' + G*Q*G';

    % save data
    Xest(i,:) = X_post; 
end
vel_est = Xest(:,2);
general_disturbance_est = Xest(:,3);

%% plot
figure;
ax = createSubplot(3,2);
hold(ax,'on')
grid(ax,'on')

plot(ax(1,1), t, pos_measured)
plot(ax(2,1), t, vel_true)
plot(ax(2,1), t, vel_est)
plot(ax(3,1), t, u_input)
plot(ax(3,2), t, d_true)
plot(ax(3,2), t, general_disturbance_est)

title(ax(1,1),'pos')
title(ax(2,1),'vel')
title(ax(3,1),'u')
title(ax(3,2),'general disturbance')
ylabel(ax(1,1),'rad')
ylabel(ax(2,1),'rad/s')
ylabel(ax(3,:),'volt')
xlabel(ax,'time (sec)')
legend(ax(1,1),'measured')
legend(ax(2,1),{'true', 'estimated'})
legend(ax(3,1),'input')
legend(ax(3,2),{'true', 'estimated'})
linkaxes(ax,'x')
loose_ylim(ax)

%% subfunction
function H = phpX(X,parm)
    % y = h(X) = x1
    % phpX = [px1px1 px1px2 px1px3];
    H = [1 0 0];
end

function F = pfpX(X,u,dt,parm)
    % X[k+1] = f(X[k],u[k]) = [f1; f2; f3]
    % pfpX = [p(f1)px1, p(f1)px2, p(f1)px3;
    %         p(f2)px1, p(f2)px2, p(f2)px3;
    %         p(f3)px1, p(f3)px2, p(f3)px3]
    % f1 = x1 + x2*dt
    % f2 = x2 + (-a*x2 + b*u - b*d_Coulomb + b*x3)*dt
    % f3 = x3
    % d_Coulomb = d_Coulomb_coeff*tanh(x2/d_Coulomb_threshold);
    % tanh(x/a)' = sech^2(x/a)/a

    % substitute
    a = parm.a;  % viscous friction [N-m / rad/s]
    b = parm.b;  % gain [N-m/volt]
    d_Coulomb_coeff = parm.d_Coulomb_coeff;         % Coulomb friction coeff [volt]
    d_Coulomb_threshold = parm.d_Coulomb_threshold; % Coulomb friction thrshould [rad/s]
    x2 = X(2);   % vel [rad/s]

    % pfpX
    F = [1, dt, 0; 
         0, 1-a*dt-b*dt*d_Coulomb_coeff/d_Coulomb_threshold*sech(x2/d_Coulomb_threshold)^2, b*dt; 
         0, 0, 1];
end

function X = plant_f(X,u,dt,parm)
    % substitute
    a = parm.a;  % viscous friction [N-m / rad/s]
    b = parm.b;  % gain [N-m/volt]
    d_Coulomb_coeff = parm.d_Coulomb_coeff;         % Coulomb friction coeff [volt]
    d_Coulomb_threshold = parm.d_Coulomb_threshold; % Coulomb friction thrshould [rad/s]
    x1 = X(1);  % pos [rad]
    x2 = X(2);  % vel [rad/s]
    x3 = X(3);  % disturbance [volt]
    
    % disturbance
    d_Coulomb = d_Coulomb_coeff*tanh(x2/d_Coulomb_threshold);
    
    % plant & integral (backward Eulaer method)
    x1_next = x1 + x2*dt;
    x2_next = x2 + (-a*x2 + b*u - b*d_Coulomb + b*x3)*dt;
    x3_next = x3;
    X = [x1_next; x2_next; x3_next];
end

function y = plant_h(X,parm)
    y = X(1);
end

function d = get_d_true(t,X,u_input,plant_parm, model_parm)
    ap = plant_parm.a;  % viscous friction [N-m / rad/s]
    bp = plant_parm.b;  % gain [N-m/volt]
    dp_Coulomb_coeff = plant_parm.d_Coulomb_coeff;         % Coulomb friction coeff [volt]
    dp_Coulomb_threshold = plant_parm.d_Coulomb_threshold; % Coulomb friction thrshould [rad/s]

    am = model_parm.a;  % viscous friction [N-m / rad/s]
    bm = model_parm.b;  % gain [N-m/volt]
    dm_Coulomb_coeff = model_parm.d_Coulomb_coeff;         % Coulomb friction coeff [volt]
    dm_Coulomb_threshold = model_parm.d_Coulomb_threshold; % Coulomb friction thrshould [rad/s]

    d = zeros(length(t),1);
    for i = 1:length(t)
        x1 = X(i,1);
        x2 = X(i,2);
        u = u_input(i);
        dp_Coulomb = dp_Coulomb_coeff*tanh(x2/dp_Coulomb_threshold);
        dm_Coulomb = dm_Coulomb_coeff*tanh(x2/dm_Coulomb_threshold);
        d(i) = ( -(ap-am)*x2 + (bp-bm)*u - ( bp*dp_Coulomb - bm*dm_Coulomb ) )/bp;
    end
end



