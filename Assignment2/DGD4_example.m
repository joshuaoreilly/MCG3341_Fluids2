
%% Parameters

Uo      = 1;                % Speed constant
U0      = 0;                % Initial speed
nu      = 10^(-6);          % Water viscosity
h       = 0.01;             % Thickness between the two plates
npts    = 20;               % Numbers of points between the two plates (mesh)
Re      = [0.1,1,100];      % Reynolds number
omega   = Re*nu/(h^2);      % Omega

%% Discretization

dy      = h/(npts-1);       % Mesh increment
dt      = 0.5*dy^2/nu;      % Maximum time increment
C  = (nu*dt)/(dy^2);   % Stability parameter C

t                = 0;
U1(npts)         = 0;  %U1, U2, U3 for Re=0.1,1,100
U2(npts)         = 0;  %Initial value of U1,U2,U3 is zero
U3(npts)         = 0;
Unew1            = U1; %Updated U1,U2,U3
Unew2            = U2;
Unew3            = U3;
y                = linspace(0,h,npts);

tsteps = 100;
tfinal = dt*tsteps;

for i = 1:tsteps

    t   = t+dt;
    Up1  = Uo*sin(omega(1)*t); % Velocity of the bottom plate, Re=0.1
    Up2  = Uo*sin(omega(2)*t); % Velocity of the bottom plate, Re=1
    Up3  = Uo*sin(omega(3)*t); % Velocity of the bottom plate, Re=100
    Unew1(1)     = Up1; % The first point: velocity=the bottom plate
    Unew1(npts)  = U0;  % The top poiint: velocity= the top plate
    Unew2(1)     = Up2;
    Unew2(npts)  = U0;
    Unew3(1)     = Up3;
    Unew3(npts)  = U0;

    for j = 2:npts-1
    
        Unew1(j) = U1(j) + C * (U1(j-1)+U1(j+1)-2*U1(j));
        Unew2(j) = U2(j) + C * (U2(j-1)+U2(j+1)-2*U2(j));
        Unew3(j) = U3(j) + C * (U3(j-1)+U3(j+1)-2*U3(j));
        j       = j+1;

    end

    figure(1)
    title(['Velocity profile for T = ',num2str(t)])
    subplot(1,3,1)
    plot(y,U1,'k','LineWidth',2)
    title(['\omega = ',num2str(omega(1))])
    xlabel('y(m)')
    ylabel('u (m/s)')
    %legend(['t = ',num2str(t),' sec'])
    axis([0,h,-Uo,Uo])
    subplot(1,3,2)
    plot(y,U2,'b','LineWidth',2)
    title(['\omega  = ',num2str(omega(2))])
    xlabel('y(m)')
    ylabel('u (m/s)')
    %legend(['t = ',num2str(t),' sec'])
    axis([0,h,-Uo,Uo])
    subplot(1,3,3)
    plot(y,U3,'r','LineWidth',2)
    title(['\omega = ',num2str(omega(3))])
    xlabel('y(m)')
    ylabel('u (m/s)')
    %legend(['t = ',num2str(t),' sec'])
    axis([0,h,-Uo,Uo])
    drawnow
    %pause(0.001);

    U1 = Unew1;
    U2 = Unew2;
    U3 = Unew3;
    %fprintf('%.f\n %',i);
end
pause(1)
