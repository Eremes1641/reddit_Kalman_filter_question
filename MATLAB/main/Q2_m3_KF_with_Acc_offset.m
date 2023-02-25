clc; clear; close all;

%% declare
Fs = 500;   % sampling freq [Hz]
dt = 1/Fs;  % sampling period [sec]

%% ODE
% simulation
time = 0:dt:1.5;
IC = [0 0]';
equ = @(t,IC)equ_plant(t,IC);
[t_sim,X_sim] = ode45(equ, time, IC);

% extract data
simData = extract_sim_data(equ, t_sim, X_sim);
t = simData.t;
pos_true = simData.x1;
vel_true = simData.x2;
acc_true = simData.acc;
u_input = simData.u;
lenT = length(t);

%% measurement data
% discrete position measurement
resolution = 1/100; % 100 CPR
pos_measured = quantize_signal(pos_true,resolution);
% add offset & noise into Acc
acc_offset = 2*ones(length(acc_true),1);
acc_measured = acc_true + acc_offset;
acc_measured = awgn(acc_measured,10,'measured');

%% kalman filter
% model
Aest = [0 1 0; 0 0 -1; 0 0 0];
Best = [0 1 0; 0 1 0; 0 0 1]';
Cest = [1 0 0];
Pest = ss(Aest,Best,Cest,[]);
Pest.StateName = {'pos', 'vel', 'acc_offset'};
Pest.InputName = {'acc', 'w_acc', 'w_acc_offset'};
Pest.OutputName = {'pos'};
Pest % display model
Pest = c2d(Pest,dt);

% kalman filter
Q = diag([1E3 1E6]);    % Acc weighting, Acc offset weighting
R = 1E0;                % pos weighting
[kalmf,~,~] = kalman(Pest,Q,R,0,1,1);
kalmf % display observer

%% run Kalman filter
obv_u = acc_measured;
obv_Y = pos_measured;
obv_IC = zeros(length(Aest),1);
[~,Xest,~] = sim_obv(kalmf,obv_u,obv_Y,obv_IC);
vel_est = Xest(:,2);
acc_offset_est = Xest(:,3);

%% plot
figure;
ax = createSubplot(3,2);
hold(ax,'on')
grid(ax,'on')

plot(ax(1,1), t, pos_measured)
plot(ax(2,1), t, vel_true)
plot(ax(2,1), t, vel_est)
plot(ax(3,1), t, acc_measured)
plot(ax(3,2), t, acc_offset)
plot(ax(3,2), t, acc_offset_est)

title(ax(1,1),'pos')
title(ax(2,1),'vel')
title(ax(3,1),'u (acc)')
title(ax(3,2),'acc offset')
ylabel(ax(1,1),'rad')
ylabel(ax(2,1),'rad/s')
ylabel(ax(3,:),'rad/s/s')
xlabel(ax,'time (sec)')
legend(ax(1,1),'measured')
legend(ax(2,1),{'true', 'estimated'})
legend(ax(3,1),'measured')
legend(ax(3,2),{'true', 'estimated'})
linkaxes(ax,'x')
loose_ylim(ax)



