clear; clc; close all;

%% Generate measured data from simulation
%initial constants/conditions to generate data
gpos = 10;
x0 = 0;
y0 = 0;
vx0 = 15;
vy0 = 12;
B = 0.9;
g = -10;
tf = roots([g/2, vy0, y0]);
tfp = find(tf > 0);
tf = tf(tfp);
xf = x0 + vx0*tf;
yf = y0 + vy0*tf + g*tf.^2/2;
max_t = 1;



%generate randomized trial x and y data from state data (with error)
numtrials = 1;
datapoints = 3;
x0data = x0 + 0.0*randn(numtrials, 1);
y0data = y0 + 0.0*randn(numtrials, 1);
vx0data = vx0 + 0.0*randn(numtrials, 1);
vy0data = vy0 + 0.0*randn(numtrials, 1);
tdata = zeros(numtrials*datapoints, 1);
xdata = zeros(numtrials*datapoints, 1);
ydata = zeros(numtrials*datapoints, 1);
for i = 1:numtrials
    %forward propagate to obtain state data at various times t
    initial_conditions = [x0data(i,:); y0data(i,:); vx0data(i,:); vy0data(i,:); B];
    [t, s] = ode45('Equations', [0, max_t], initial_conditions);
    
    for j = 1:datapoints
        curr_index = j+(i-1)*datapoints;
        tdata(curr_index,:) = t(round(size(t,1)*(j/datapoints)),:);
        xdata(curr_index,:) = s(round(size(t,1)*(j/datapoints)),1);
        ydata(curr_index,:) = s(round(size(t,1)*(j/datapoints)),2);
    end
end

mdata = zeros(numtrials*datapoints*2, 1);
mdata(1:2:end) = xdata;
mdata(2:2:end) = ydata;

figure(1);
[t, traj] = ode45('Equations', [0, max_t], initial_conditions);
plot(t, traj(:,2), 'g');
hold on;

%% Estimation starts here
dB = 0.001;

% H = zeros(numtrials*datapoints*2, 4);
% H(1:2:end, 1) = 1;
% H(2:2:end, 2) = 1;
% H(1:2:end, 3) = tdata;
% H(2:2:end, 4) = tdata;
% ym = zeros(numtrials*datapoints*2, 1);
% ym(1:2:end) = xdata;
% ym(2:2:end) = ydata + g*(tdata.^2)/2;
% PinH = pinv(H);
% xest = PinH*ym;

%Initial Guess:
x0_e = 0;
y0_e = 0;
vx0_e = 15;
vy0_e = 12;
B_e = 2.0;
delta_B = 1;
toler = 0.001;

initial_conditions_e = [x0_e; y0_e; vx0_e; vy0_e; B_e];
[t, traj] = ode45('Equations', [0, max_t], initial_conditions_e);
figure(1)
plot(t, traj(:,2));
hold on;

convergence = zeros(500, 1);
convergence(1) = B_e;
i = 0;

deltas = zeros(500, 1);

while abs(delta_B) > toler && abs(delta_B) < 10
    
    initial_conditions_e = [x0_e; y0_e; vx0_e; vy0_e; B_e];
    [t_est, s_est] = ode45('Equations', [0, max_t], initial_conditions_e);

    xdata_est = interp1(t_est, s_est(:,1), tdata);
    ydata_est = interp1(t_est, s_est(:,2), tdata);

    estdata = zeros(numtrials*datapoints*2, 1);
    estdata(1:2:end) = xdata_est;
    estdata(2:2:end) = ydata_est;

    error1 = (mdata - estdata).^2;
    cost1 = sum(error1);
    
    %used to find change in cost over small nudge to B
    initial_conditions_e = [x0_e; y0_e; vx0_e; vy0_e; B_e+dB];
    [t_est, s_est] = ode45('Equations', [0, max_t], initial_conditions_e);

    xdata_est = interp1(t_est, s_est(:,1), tdata);
    ydata_est = interp1(t_est, s_est(:,2), tdata);

    estdata = zeros(numtrials*datapoints*2, 1);
    estdata(1:2:end) = xdata_est;
    estdata(2:2:end) = ydata_est;

    error2 = (mdata - estdata).^2;
    cost2 = sum(error2);

    approx_derivative = (cost2-cost1)/(dB);

    %delta_B = -(cost1)./approx_derivative;
    %B_e = B_e + delta_B;
    delta_B = -(cost1)./approx_derivative;
    B_e = B_e + delta_B;
    i = i + 1;
    convergence(i+1) = B_e;
    deltas(i+1) = delta_B;
end

[t, traj] = ode45('Equations', [0, max_t], initial_conditions_e);
plot(t, traj(:,2));
legend("Initial guess", "After convergence");
title("Projectile trajectories of intial guess and after convergence");
% if B_e < 0
%     B_e = B_e * -1;
% end

convergence = nonzeros(convergence);
figure(2);
plot(0:i,convergence);
title("Convergence of estimated parameter");

deltas = nonzeros(deltas);
figure(3);
plot(0:i-1,deltas);
title("Convergence of deltas");





