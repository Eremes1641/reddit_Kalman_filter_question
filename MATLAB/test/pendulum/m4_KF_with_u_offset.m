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
model_parm.a = 0.5;  % viscous friction [N-m / rad/s]
model_parm.b = 5;    % gain [N-m/volt]
model_parm.g = 9.81; % gravity
model_parm.l = 0.1;  % length [m]

% Kalman filter
Q = diag([1E3 1E10]); % u weighting, u offset weighting
R = 1E0;              % pos weighting

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
d_true = simData.d;
u_input = simData.u;
lenT = length(t);

%% measurement data
% discrete position measurement
resolution = 1/100; % 100 CPR
pos_measured = quantize_signal(pos_true,resolution);

%% kalman filter
% model
a = model_parm.a;  % viscous friction [N-m / rad/s]
b = model_parm.b;  % gain [N-m/volt]
Aest = [0 1 0; 0 -a -1; 0 0 0];
Best = [0 b 0; 0 1 0; 0 0 1]';
Cest = [1 0 0];
Pest = ss(Aest,Best,Cest,[]);
Pest.StateName = {'pos', 'vel', 'd'};
Pest.InputName = {'u', 'w_u', 'w_d'};
Pest.OutputName = {'pos'};
Pest % display model
Pest = c2d(Pest,dt);

% kalman filter
[kalmf,~,~] = kalman(Pest,Q,R,0,1,1);
kalmf % display observer

%% run Kalman filter
obv_u = u_input;
obv_Y = pos_measured;
obv_IC = zeros(length(Aest),1);
[~,Xest,~] = sim_obv(kalmf,obv_u,obv_Y,obv_IC);
vel_est = Xest(:,2);
u_offset_est = Xest(:,3);

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
plot(ax(3,2), t, u_offset_est)

title(ax(1,1),'pos')
title(ax(2,1),'vel')
title(ax(3,1),'u')
title(ax(3,2),'general disturbance')
ylabel(ax(1,1),'rad')
ylabel(ax(2,1),'rad/s')
ylabel(ax(3,:),'N-m')
xlabel(ax,'time (sec)')
legend(ax(1,1),'measured')
legend(ax(2,1),{'true', 'estimated'})
legend(ax(3,1),'input')
legend(ax(3,2),{'true', 'estimated'})
linkaxes(ax,'x')
loose_ylim(ax)


