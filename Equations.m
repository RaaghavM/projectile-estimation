%%%
%% This file contains the system dynamics

function dots = dynamics(t, s)

%%%
%% Define system parameters
g = 9.8; % gravitational acceleration constant
p = 1.2; % flow density which is approx atmospheric density at low altitudes
%m = 1; % mass of projectile
%S = 0; %0.25 %0.8;
%C_D = 0.01;
%C_L = 0; %0.01; %0.5; % lift coefficient
%B = 0.9; %m/(S*C_D); % ballistic coefficient - estimates ratio of gravity to drag
%%%
%% Provide time derivatives of the states. ODE45 will use these to
%% perform the integration
x = s(1);
y = s(2);
v_x = s(3);
v_y = s(4);
B = s(5);

xdot =  v_x;
ydot = v_y;

v = sqrt(v_x^2 + v_y^2);

v_xdot = -1*(p/(2*B))*v*v_x; %-(1/(2*m))*p*v*S*C_L*v_y;
v_ydot = -1*(g+(p/(2*B))*v*v_y); %+(1/(2*m))*p*v*S*C_L*v_x;

%%%
%% collect all derivatives
%%%
dots = [xdot; ydot; v_xdot; v_ydot; 0];    %This must also be a COLUMN vector.