clc; clear all; close all

% Setup

U0 = 0;     % Initial speed
nu = 10^(-6);           % Kinematic viscocity nu = mu/rho
h = 100;               % Height
nmbrPts = 20;           % Number of vertical divides
% Reynolds Number
Re = 0.01;        % Value of Re changed for each question
omega   = Re*nu/(h^2);      % Omega

% Discretization
dy = h/(nmbrPts - 1);      % height step --- MINUS ONE FOR SOME REASON
dt = 0.5 * (dy^2) / nu; % time step
C = (nu*dt)/(dy^2);     % Stability parameter

elapsedTime = int32(10/omega);        % Measure profile at this time
                                        % Re 0.001 was chosen because it requires much more time to become stable
tsteps = idivide(elapsedTime, dt, 'ceil');  % Find number of steps, round up

t = 0;
Uo = 1;     % Speed constant
U(nmbrPts) = 0;    % Velocities for first Reynolds number
% Rows are the speeds for different times, columns are at different heights
%Unew = zeros(nmbrPts + 2);  % Updating velocity values or something...
y = linspace(0,h,nmbrPts);

for i = 1:tsteps
  %Uplate = Uo*sin(omega(1)*t);   % Velocity of bottom plate
  t = t+dt;
  for j = 2:nmbrPts-1
    U(j) = U(j) + C * (U(j-1) + U(j+1) - (2*U(j)));  % Update speeds for current time - First Re
  end
  Uplate = cos(omega*t);      % Set plate speed
  U(end) = U0;                % Set top plate speed
end

plot(U, y, 'k');
title('Velocity profile for R = 0.01');
xlabel('Velocity (m/s)');
ylabel('height (m)')
grid minor
