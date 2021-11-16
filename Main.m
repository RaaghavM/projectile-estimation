clear all; clc; %%close all;
%% Specify initial conditions
%%%
theta = pi/4;
v0mag = 30;
if theta ~= -1 && v0mag ~= -1
   vx0 = v0mag*cos(theta);
   vy0 = v0mag*sin(theta);
else
    vx0 = 5;
    vy0 = 5;
end

initial_conditions = [0; 0; vx0; vy0]; %initial posx, posy, vx, vy
max_t = 10
                                
%%%
%% Call ode45 to integrate the system of equations
%%%
[t s] = ode45('Equations', [0, max_t], initial_conditions);


%%%
%% Plot the time histories of the states
%%%
%figure(1)
%plot(t, s(:,1), 'k', 'linewidth', 2);
%hold on;
%plot(t, s(:,2), 'k-.', 'linewidth', 2);
%set(gca, 'fontsize', 12, 'fontweight', 'bold');
%xlabel('Time, s');
%ylabel('State');
%legend('x(t)', 'v(t)');
%grid on;

figure(1)
plot(s(:,1), s(:,2), 'r', 'linewidth', 2);
hold on;
plot(s(:,3), s(:,4), 'r.', 'linewidth', 2);
set(gca, 'fontsize', 12, 'fontweight', 'bold');
xlabel('x');
ylabel('y');
legend({'Position', 'Velocity'});
title('Trajectory');
ylim([0,max(s(:,2))])
xlim([0,max(s(:,1))])
grid on;