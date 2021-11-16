clear all; clc; close all;

g = -9.8;
dy = -1.6;

D1 = readmatrix('Projectile Experiment Data - Day 1.csv')
D2 = readmatrix('Projectile Experiment Data - Day 2.csv')
D3 = readmatrix('Projectile Experiment Data - Day 3.csv')
D4 = readmatrix('Projectile Experiment Data - Day 4.csv') %%lessened force of rubber band
D5 = readmatrix('Projectile Experiment Data - Day 5.csv') %%ping pong ball (lessened force?)
PerfectTheoretical = readmatrix('Projectile Experiment Data - Theoretical.csv')
RandomTheoretical = readmatrix('Projectile Experiment Data - RandomTheoretical.csv')
D_All = [D1;D2;D3]

v_x0 = D1(:,4) ./ D1(:,2);
D1 = [D1 v_x0];
v_x0 = D2(:,4) ./ D2(:,2);
D2 = [D2 v_x0];
v_x0 = D3(:,4) ./ D3(:,2);
D3 = [D3 v_x0];
v_x0 = D4(:,4) ./ D4(:,2);
D4 = [D4 v_x0];
v_x0 = D5(:,4) ./ D5(:,2);
D5 = [D5 v_x0];
v_x0 = D_All(:,4) ./ D_All(:,2);
D_All = [D_All v_x0];


v_y0 = (dy - (0.5 * g) .* (D1(:,2) .* D1(:,2))) ./ D1(:,2);
D1 = [D1 v_y0];
v_y0 = (dy - (0.5 * g) .* (D2(:,2) .* D2(:,2))) ./ D2(:,2);
D2 = [D2 v_y0];
v_y0 = (dy - (0.5 * g) .* (D3(:,2) .* D3(:,2))) ./ D3(:,2);
D3 = [D3 v_y0];
v_y0 = (dy - (0.5 * g) .* (D4(:,2) .* D4(:,2))) ./ D4(:,2);
D4 = [D4 v_y0];
v_y0 = (dy - (0.5 * g) .* (D5(:,2) .* D5(:,2))) ./ D5(:,2);
D5 = [D5 v_y0];
v_y0 = (dy - (0.5 * g) .* (D_All(:,2) .* D_All(:,2))) ./ D_All(:,2);
D_All = [D_All v_y0];

%v_x0 = zeros(size(D1,1),1);
%h = histogram(M(:,2),10)
figure(1)
h1 = histogram(D5(:,4), 'normalization', 'probability')
figure(2)
h2 = histogram(D1(:,4), 'normalization', 'probability')
%scatter(D_All(:,),D_All(:,6))

numdata = size(D1, 1);

gpos = 9.81;
tf = D1(:,2);
xf = D1(:,4);
yf = 0;
%H = zeros(2*numdata, 4);
%H = zeros(2*numdata, 3);
H = zeros(2*numdata, 4);
H(1:2:end, 1) = 1;
H(2:2:end, 2) = 1;
H(1:2:end, 3) = tf;
H(2:2:end, 4) = tf;

ym = zeros(2*numdata, 1);
ym(1:2:end) = xf;
ym(2:2:end) = yf + 0.5*gpos* (tf .* tf); %should multiply tf by its mean?

PinH = pinv(H);
xest = PinH*ym

figure(3)
initial_conditions = xest;
[t, traj] = ode45('Equations', [0, 2], initial_conditions);
plot(t(:,1),  traj(:,2));

figure(4)
err = ym - H*xest;
plot(1:2*numdata, err, 'bo');
sum(err .* err);

mean(D1(:,2)); %mean tf
mean(D1(:,4)); %mean xf

% gpos = 9.81;
% tf = D1(:,2);
% xf = D1(:,4);
% yf = 0;
% H = zeros(2*numdata, 3);
% H(1:2:end, 1) = 1;
% H(2:2:end, 2) = tf;
% H(1:2:end, 3) = tf;
% 
% ym = zeros(2*numdata, 1);
% ym(1:2:end) = xf;
% ym(2:2:end) = yf + 0.5*gpos* (tf .^ 2) + dy; 
% 
% PinH = pinv(H);
% xest = PinH*ym
% 
% figure(3)
% initial_conditions = xest;
% [t, traj] = ode45('Equations', [0, 2], initial_conditions);
% plot(t(:,1),  traj(:,2));
% 
% figure(4)
% err = ym - H*xest;
% plot(1:2*numdata, err, 'bo');
% sum(err .* err)
% 
% mean(D1(:,2)) %mean tf
% mean(D1(:,4)) %mean xf

% gpos = 10
% x0 = 0;
% y0 = 0;
% vx0 = 10;
% vy0 = 6;
% g = -10;
% tf = roots([g/2, vy0, y0]);
% tfp = find(tf > 0);
% tf = tf(tfp);
% 
% xf = x0 + vx0*tf;
% yf = y0 + vy0*tf + g*tf^2/2;
% 
% numdata = 50;
% xfdata_unn = randn(numdata, 1);
% tfdata_unn = randn(numdata, 1);
% 
% xfdata = xf + 1*xfdata_unn;
% tfdata = tf + 0.1*tfdata_unn;
% 
% yf = 0;
% H = zeros(2*numdata, 4);
% %H = zeros(2*numdata, 3);
% %H = zeros(2*numdata, 4);
% H(1:2:end, 1) = 1;
% H(2:2:end, 2) = 1;
% H(1:2:end, 3) = tfdata;
% H(2:2:end, 4) = tfdata;
% 
% ym = zeros(2*numdata, 1);
% ym(1:2:end) = xfdata;
% ym(2:2:end) = yf + 0.5*gpos* (tfdata .^ 2); 
% 
% PinH = pinv(H);
% xest = PinH*ym
% 
% figure(3)
% initial_conditions = xest;
% [t, traj] = ode45('Equations', [0, 2], initial_conditions);
% plot(t(:,1),  traj(:,2));
% 
% figure(4)
% err = ym - H*xest;
% plot(1:2*numdata, err, 'bo');
% sum(err .* err)
% 
% mean(D1(:,2)) %mean tf
% mean(D1(:,4)) %mean xf


