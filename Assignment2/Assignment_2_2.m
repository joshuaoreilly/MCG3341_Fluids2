clc; clear all; close all;

%% Parameters

Uo      = 1;                % Speed constant
U0      = 0;                % Initial speed
nu      = 10^(-6);          % Water viscosity
h       = 0.01;             % Thickness between the two plates
npts    = 20;               % Numbers of points between the two plates (mesh)
Re      = [0.01,100];      % Reynolds number
omega   = Re*nu/(h^2);      % Omega

%% Discretization

dy      = h/(npts-1);       % Mesh increment
dt      = 0.5*dy^2/nu;      % Maximum time increment
C  = (nu*dt)/(dy^2);   % Stability parameter C

t                = 0;
U1(npts)         = 0;
U2(npts)         = 0;
Unew1            = U1;
Unew2            = U2;
y                = linspace(0,h,npts);

%tsteps = 10/omega(1);
tsteps = 100000;         %took about 140s => t=10000
                         %ran for about 6 minutes with t=100000 and no
                         %result
tsteps1= 10/omega(2);
tfinal = dt*tsteps;

for i = 1:tsteps

    t   = t+dt;
    Up1  = Uo*cos(omega(2)*t);
    Unew1(1)     = Up1;
    Unew1(npts)  = U0;

    for j = 2:npts-1

        Unew1(j) = U1(j) + C * (U1(j-1)+U1(j+1)-2*U1(j));

        j       = j+1;

    end
    %fprintf('%.f\n %',i);
end

figure(1)
title(['Velocity profile for T = ',num2str(t)])
subplot(1,2,1)
plot(y, U1,'k','LineWidth',2)
title(['\omega = ',num2str(omega(1))])
xlabel('y(m)')
ylabel('u (m/s)')
grid

axis([0,h,-Uo,Uo])

U1 = Unew1;

%for i = 1:tsteps1

    %t   = t+dt;
    %Up2  = Uo*cos(omega(2)*t);
    %Unew2(1)     = Up2;
    %Unew2(npts)  = U0;

    %for j = 2:npts-1

        %Unew2(j) = U2(j) + C * (U2(j-1)+U2(j+1)-2*U2(j));

        %j       = j+1;

    %end

    %axis([0,h,-Uo,Uo])
    %subplot(1,2,2)
    %plot(y,U2,'b','LineWidth',2)
    %title(['\omega  = ',num2str(omega(2))])
    %xlabel('y(m)')
    %ylabel('u (m/s)')
    %grid
    %legend(['t = ',num2str(t),' sec'])
    %axis([0,h,-Uo,Uo])
    %pause(0.001);

    %U2 = Unew2;
    %fprintf('%.f\n %',i);
%end
