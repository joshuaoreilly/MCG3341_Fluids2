clc; clear all; close all

% Setup

U0 = 0;     % Initial speed
nu = 10^(-6);           % Kinematic viscocity nu = mu/rho
h = 0.01;               % Height
nmbrPts = 20;           % Number of vertical divides
% Reynolds Number?
Re = 0.01;
omega   = Re*nu/(h^2);      % Omega

% Discretization
dy = h/(nmbrPts - 1);      % height step --- MINUS ONE FOR SOME REASON
dt = 0.5 * (dy^2) / nu; % time step
C = (nu*dt)/(dy^2);     % Stability parameter

elapsedTime = int32(10/omega);        % Measure profile at this time
tsteps = idivide(elapsedTime, dt, 'ceil');  % Find number of steps, round up

t = 0;
Uo = 1;     % Speed constant
U = zeros(tsteps, nmbrPts);     % Velocities
% Rows are the speeds for different times, columns are at different heights
%Unew = zeros(nmbrPts + 2);  % Updating velocity values or something...
y = linspace(0,h,nmbrPts);

% FIRST ITERATION
Uplate = cos(omega*t);
Uup = [U(1,2:end) 0];
Udown = [Uplate U(1,1:end-1)];
U(1, :) = U(1, :) + (C .* (Udown + Uup - 2.*U(1, :)));


for i = 2:tsteps
  %Uplate = Uo*sin(omega(1)*t);   % Velocity of bottom plate
  Uplate = cos(omega*t);
  % [1 2 3 4] becomes [2 3 4 0]
  Uup = [U((i-1),2:end) 0];
  % [1 2 3 4] becomes [Uplate 1 2 3]
  Udown = [Uplate U(i-1,1:end-1)];
  t = t+dt;
  U(i, :) = U(i-1, :) + (C .* (Udown + Uup - 2.*U(i-1, :)));  % Update speeds for current time
  U(i, 1) = Uplate;   % Set plate speed
end

plot(U(end,:), y, 'k');
title(['Velocity profile for Re = ',num2str(Re)]);
xlabel('Velocity (m/s)');
ylabel('height (m)')
grid minor
