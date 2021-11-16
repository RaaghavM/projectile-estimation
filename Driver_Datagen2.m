clear; clc; close all;

gpos = 10;
x0 = 0;
y0 = 0;
vx0 = 10;
vy0 = 6;
g = -10;
tf = roots([g/2, vy0, y0]);
tfp = find(tf > 0);
tf = tf(tfp);
xf = x0 + vx0*tf;
yf = y0 + vy0*tf + g*tf.^2/2;

%generating x data (discarding corresponding y data)
numtrials_X = 10000;
datapoints_X = 1;
x0data = x0 + 0.1*randn(numtrials_X, 1);
y0data = y0 + 0.1*randn(numtrials_X, 1);
vx0data = vx0 + 0.05*randn(numtrials_X, 1);
vy0data = vy0 + 0.05*randn(numtrials_X, 1);
tdata_X = zeros(numtrials_X*datapoints_X, 1);
xdata = zeros(numtrials_X*datapoints_X, 1);
for i = 1:numtrials_X
    tf_i = max(roots([g/2, vy0data(i), y0data(i)]));
    for j = 1:datapoints_X
        curr_index = j+(i-1)*datapoints_X;
        tdata_X(curr_index) = (j/datapoints_X)*tf_i;
        xdata(curr_index) = x0data(i) + vx0data(i)*tdata_X(curr_index);
    end
end

%estimating initial x-velocity and x-position using tf and final x-position
H_X = zeros(numtrials_X*datapoints_X, 2);
H_X(:, 1) = 1;
H_X(:, 2) = tdata_X;

ym_X = xdata;

PinH_X = pinv(H_X);
xest_X = PinH_X*ym_X


%generating y data (discarding corresponding x data)
numtrials_Y = 10000;
datapoints_Y = 2;
x0data = x0 + 0.1*randn(numtrials_Y, 1);
y0data = y0 + 0.1*randn(numtrials_Y, 1);
vx0data = vx0 + 0.05*randn(numtrials_Y, 1);
vy0data = vy0 + 0.05*randn(numtrials_Y, 1);
tdata_Y = zeros(numtrials_Y*datapoints_Y, 1);
ydata = zeros(numtrials_Y*datapoints_Y, 1);
for i = 1:numtrials_Y
   tf_i = max(roots([g/2, vy0data(i), y0data(i)]));
   for j = 1:datapoints_Y
        curr_index = j+(i-1)*datapoints_Y;
        tdata_Y(curr_index) = (j/datapoints_Y)*tf_i;
        ydata(curr_index) = y0data(i) + vy0data(i)*tdata_Y(curr_index) + g*tdata_Y(curr_index).^2/2;
   end
end

%estimating initial y-velocity and y-position using tf and final y-position
H_Y = zeros(numtrials_Y*datapoints_Y, 2);
H_Y(:, 1) = 1;
H_Y(:, 2) = tdata_Y;

ym_Y = ydata - 0.5*g*(tdata_Y.^2);

PinH_Y = pinv(H_Y);
xest_Y = PinH_Y*ym_Y

A = zeros(2, 2);
C = H_Y;
Ob = obsv(A,C);
unob = length(A)-rank(Ob)

