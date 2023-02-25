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
u_input = simData.u;
lenT = length(t);

%% measurement data
% quantize position measurement
resolution = 1/100; % 100 CPR
pos_measured = quantize_signal(pos_true,resolution);

%% difference & low pass filter
bandwidth = 100; % rad/s
vel_est = LPF(pos_measured,bandwidth,dt);

%% plot
figure;
ax = createSubplot(3,1);
hold(ax,'on')
grid(ax,'on')

plot(ax(1), t, pos_measured)
plot(ax(2), t, vel_true)
plot(ax(2), t, vel_est)
plot(ax(3), t, u_input)

title(ax(1),'pos')
title(ax(2),'vel')
title(ax(3),'u (input)')
ylabel(ax(1),'rad')
ylabel(ax(2),'rad/s')
ylabel(ax(3),'volt')
xlabel(ax,'time (sec)')
legend(ax(1),'measured')
legend(ax(2),{'true', 'estimated'})
legend(ax(3),'input')
linkaxes(ax,'x')
loose_ylim(ax)


