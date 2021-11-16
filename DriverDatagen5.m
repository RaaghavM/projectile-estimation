clear; clc; close all;

%% Generate measured data from simulation
%initial constants/conditions to generate data
gpos = 10;
x0 = 1;
y0 = 1;
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
numtrials = 100;
datapoints = 5;
uniform_datapoints = true;
x0data = x0 + 0.05*randn(numtrials, 1); %0.01 standard deviation is max error for good convergence
y0data = y0 + 0.05*randn(numtrials, 1);
vx0data = vx0 + 0.05*randn(numtrials, 1);
vy0data = vy0 + 0.05*randn(numtrials, 1);
tdata = zeros(numtrials*datapoints, 1);
xdata = zeros(numtrials*datapoints, 1);
ydata = zeros(numtrials*datapoints, 1);
for i = 1:numtrials
    %forward propagate to obtain state data at various times t
    initial_conditions = [x0data(i,:); y0data(i,:); vx0data(i,:); vy0data(i,:); B];
    [t, s] = ode45('Equations', [0, max_t], initial_conditions);
    
    for j = 1:datapoints
        if uniform_datapoints == true
            curr_index = j+(i-1)*datapoints;
            tdata(curr_index,:) = t(round(size(t,1)*(j/datapoints)),:);
            xdata(curr_index,:) = s(round(size(t,1)*(j/datapoints)),1) + 0.05*randn(1, 1);
            ydata(curr_index,:) = s(round(size(t,1)*(j/datapoints)),2) + 0.05*randn(1, 1);
        else
            curr_index = j+(i-1)*datapoints;
            tdata(curr_index,:) = t(round(size(t,1)*(j/datapoints)),:);
            xdata(curr_index,:) = s(round(size(t,1)*(j/datapoints)),1) + 0.0*randn(1, 1);
            ydata(curr_index,:) = s(round(size(t,1)*(j/datapoints)),2) + 0.0*randn(1, 1);
        end
    end
end

mdata = zeros(numtrials*datapoints*2, 1);
mdata(1:2:end) = xdata;
mdata(2:2:end) = ydata;



%% Estimation starts here

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
x0_e = 5;
y0_e = -5;
vx0_e = 13;
vy0_e = 15;
B_e = 1.3;
delta_p = [1,1,0,0,0];

toler = 0.0000001;
maxiterations = 50;

initial_conditions_e = [x0_e; y0_e; vx0_e; vy0_e; B_e];
[t, traj] = ode45('Equations', [0, max_t], initial_conditions_e);
figure(1)
plot(t, traj(:,2));
hold on;

convergence = zeros(5000, 5);
convergence(1,:) = [x0_e, y0_e, vx0_e, vy0_e, B_e];
i = 0;

deltas = zeros(5000, 5);

while (any(abs(delta_p) > toler) && all(abs(delta_p) < 100)) && (i <= maxiterations)
    
    %[t_est_c, s_est_c] = ode45('Equations', [0 max_t], initial_conditions_e);
    %xdata_est = interp1(t_est_c, s_est_c(:,1), tdata);
    %ydata_est = interp1(t_est_c, s_est_c(:,2), tdata);
    
    sol = ode45(@Equations, [0, max_t], initial_conditions_e);
    all_data = deval(sol, tdata);
    xdata_est = transpose(all_data(1,:));
    ydata_est = transpose(all_data(2,:));
    
    

    estdata1 = zeros(size(mdata));
    estdata1(1:2:end) = xdata_est;
    estdata1(2:2:end) = ydata_est;
    
    numeric_pderivs = zeros(2*numtrials*datapoints, 5);
    for j = 1:5
        numeric_pderivs(:,j) = numeric_pderiv(estdata1,initial_conditions_e,j,tdata,mdata,max_t,'Equations');
    end
    
    gradient = numeric_pderiv1(numeric_pderivs,estdata1,mdata);
    
    second_gradient = numeric_pderiv2(numeric_pderivs,numtrials,datapoints);
    
    
    delta_p = -1*(pinv(second_gradient) * gradient);
    
    
    
    initial_conditions_e = initial_conditions_e + delta_p;
    
    i = i + 1;
    convergence(i+1,:) = initial_conditions_e;
    deltas(i+1,:) = delta_p;
end

convergence(i+1,:) = [x0, y0, vx0, vy0, B];

[t, traj] = ode45('Equations', [0, max_t], initial_conditions_e);
plot(t, traj(:,2));
hold on;
initial_conditions_original = [mean(x0data); mean(y0data); mean(vx0data); mean(vy0data); B];
[to, trajo] = ode45('Equations', [0, max_t], initial_conditions_original);
plot(to, trajo(:,2));
legend("Initial guess", "After convergence", "Observed trajectory");
title("Projectile trajectories of intial guess, after convergence, and observed");

convergence = convergence(any(convergence,2),:); %remove all rows with all zeros (unused)
figure(2);
plot(0:i,convergence(:,1));
hold on;
plot(0:i,convergence(:,2));
hold on;
plot(0:i,convergence(:,3));
hold on;
plot(0:i,convergence(:,4));
hold on;
plot(0:i,convergence(:,5));
title("Convergence of estimated parameters");
legend("x0", "y0", "vx0", "vy0", "B");

deltas = deltas(any(deltas,2),:);
figure(3);
plot(0:i-1,deltas(:,1));
hold on;
plot(0:i-1,deltas(:,2));
hold on;
plot(0:i-1,deltas(:,3));
hold on;
plot(0:i-1,deltas(:,4));
hold on;
plot(0:i-1,deltas(:,5));
title("Convergence of deltas");
legend("x0", "y0", "vx0", "vy0", "B");







