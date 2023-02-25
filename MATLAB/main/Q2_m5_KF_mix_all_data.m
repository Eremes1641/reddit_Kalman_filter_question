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
d_Coulomb = simData.d;
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
a = 3;  % viscous friction [N-m / rad/s]
b = 30; % gain [N-m/volt]
Aest = [0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
Best = [0 0 0; 0 0 0; diag([1 1 1])];
Cest = [1 0 0 0 0; 0 0 1 1 0; 0 a/b 1/b 0 1];
Pest = ss(Aest,Best,Cest,[]);
Pest.StateName = {'pos', 'vel', 'acc', 'acc_offset', 'u_offset'};
Pest.InputName = {'w_acc', 'w_acc_offset', 'w_u_offset'};
Pest.OutputName = {'pos', 'acc', 'u'};
Pest % display model
Pest = c2d(Pest,dt);

% kalman filter
Q = diag([1E6 1E0 1E6]);    % acc, acc offset, u offset
R = diag([1E0 1E3 1E3]);    % pos, acc, u 
[kalmf,~,~] = kalman(Pest,Q,R,0,[1 2 3],[]);
kalmf % display observer

%% run Kalman filter
u = [];
Y = [pos_measured acc_measured u_input];
IC = zeros(length(Aest),1);
[~,Xest,~] = sim_obv(kalmf,u,Y,IC);
vel_est = Xest(:,2);
acc_est = Xest(:,3);
acc_offset_est = Xest(:,4);
u_offset_est = Xest(:,5);

%% plot
figure;
ax = createSubplot(4,2);
hold(ax,'on')
grid(ax,'on')

plot(ax(1,1), t, pos_measured)
plot(ax(2,1), t, vel_true)
plot(ax(2,1), t, vel_est)
plot(ax(3,1), t, acc_measured)
plot(ax(3,1), t, acc_est)
plot(ax(3,2), t, acc_offset)
plot(ax(3,2), t, acc_offset_est)
plot(ax(4,1), t, u_input)
plot(ax(4,2), t, d_Coulomb)
plot(ax(4,2), t, u_offset_est)

title(ax(1,1),'pos')
title(ax(2,1),'vel')
title(ax(3,1),'acc')
title(ax(3,2),'acc offset')
title(ax(4,1),'u')
title(ax(4,2),'u offset')
ylabel(ax(1,1),'rad')
ylabel(ax(2,1),'rad/s')
ylabel(ax(3,:),'rad/s/s')
ylabel(ax(4,:),'volt')
xlabel(ax,'time (sec)')
legend(ax(1,1),'measured')
legend(ax(2,1),{'true', 'estimated'})
legend(ax(3,1),{'measured', 'estimated'})
legend(ax(3,2),{'true', 'estimated'})
legend(ax(4,1),'input')
legend(ax(4,2),{'true', 'estimated'})
linkaxes(ax,'x')
loose_ylim(ax)


