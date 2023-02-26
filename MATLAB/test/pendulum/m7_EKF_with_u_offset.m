clc; clear; close all;

%% declare
% simulation
Fs = 500;   % sampling freq [Hz]
dt = 1/Fs;  % sampling period [sec]

% plant parameter
plant_parm.a = 0.5;  % viscous friction [N-m / rad/s]
plant_parm.b = 5;    % gain [N-m/volt]
plant_parm.g = 9.81; % gravity
plant_parm.l = 0.1;  % length [m]

% model parameter
model_parm.a = 0.4;  % viscous friction [N-m / rad/s]
model_parm.b = 6;    % gain [N-m/volt]
model_parm.g = 9.81; % gravity
model_parm.l = 0.15; % length [m]

% EKF parameter
Q = diag([1E3 1E7]); % u weighting, general disturbance weighting
R = 1E3;             % pos weighting

%% ODE
% simulation
time = 0:dt:1.5;
IC = [0 0]';
equ = @(t,IC)plant_pendulum(t,IC,plant_parm);
[t_sim,X_sim] = ode45(equ, time, IC);

% extract data
simData = extract_sim_data(equ, t_sim, X_sim);
t = simData.t;
pos_true = simData.x1;
vel_true = simData.x2;
u_input = simData.u;
lenT = length(t);
% d_true = simData.d;
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
ylabel(ax(3,1),'volt')
ylabel(ax(3,2),'N-m')
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
    % f2 = x2 + (-a*x2 + b*u - g/l*sin(x1) + x3)*dt
    % f3 = x3

    % substitute
    a = parm.a; % viscous friction [N-m / rad/s]
    g = parm.g; % gravity
    l = parm.l; % length [m]
    x1 = X(1);  % pos [rad]

    % pfpX
    F = [1, dt, 0;
        -g/l*cos(x1)*dt, 1-a*dt, dt;
        0, 0, 1];
end

function X = plant_f(X,u,dt,parm)
    % substitute
    a = parm.a; % viscous friction [N-m / rad/s]
    b = parm.b; % gain [N-m/volt]
    g = parm.g; % gravity
    l = parm.l; % length [m]
    x1 = X(1);  % pos [rad]
    x2 = X(2);  % vel [rad/s]
    x3 = X(3);  % disturbance [volt]
    
    % disturbance
    d_gravity = g/l*sin(x1);
    
    % plant & integral (backward Eulaer method)
    x1_next = x1 + x2*dt;
    x2_next = x2 + (-a*x2 + b*u - d_gravity + x3)*dt;
    x3_next = x3;
    X = [x1_next; x2_next; x3_next];
end

function y = plant_h(X,parm)
    y = X(1);
end

function d = get_d_true(t,X,u_input,plant_parm, model_parm)
    ap = plant_parm.a;
    bp = plant_parm.b;
    gp = plant_parm.g;
    lp = plant_parm.l;

    am = model_parm.a;
    bm = model_parm.b;
    gm = model_parm.g;
    lm = model_parm.l;
    
    d = zeros(length(t),1);
    for i = 1:length(t)
        x1 = X(i,1);
        x2 = X(i,2);
        u = u_input(i);
        d(i) = -(ap-am)*x2 + (bp-bm)*u - ( gp/lp*sin(x1) - gm/lm*sin(x1));
    end
end



