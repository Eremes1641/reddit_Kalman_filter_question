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
d_Coulomb = simData.d;
u_input = simData.u;
lenT = length(t);

%% measurement data
% discrete position measurement
resolution = 1/100; % 100 CPR
pos_measured = quantize_signal(pos_true,resolution);

%% kalman filter
% model
a = model_parm.a;
b = model_parm.b;
Aest = [0 1 0; 0 -a b; 0 0 0];
Best = [0 1 0; 0 1 0; 0 0 1]';
Cest = [1 0 0];
Pest = ss(Aest,Best,Cest,[]);
Pest.StateName = {'pos', 'vel', 'u_offset'};
Pest.InputName = {'u', 'w_u', 'w_u_offset'};
Pest.OutputName = {'pos'};
Pest % display model
Pest = c2d(Pest,dt);

% kalman filter
Q = diag([1E3 1E6]);    % u weighting, u offset weighting
R = 1E0;                % pos weighting
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
plot(ax(3,2), t, d_Coulomb)
plot(ax(3,2), t, u_offset_est)

title(ax(1,1),'pos')
title(ax(2,1),'vel')
title(ax(3,1),'u')
title(ax(3,2),'u offset')
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


