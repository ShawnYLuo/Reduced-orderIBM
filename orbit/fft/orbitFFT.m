% In this example, we look at an elliptical membrane in a fluid that is
% initially at rest.  The ellipse should oscillate and eventually settle to
% a circle with a prescribed radius.

clear;
% Add PATH reference in order to run solver
addpath('./utils');

% The number of grid points.
Nx = 4*round(2.^4);
Ny = Nx;

% Parameter values.
mu = 1;        % Viscosity.
sigma = 4e3;     % Bending constant.
rho = 1;       % Density.
rmin = 0.2;    % Length of semi-minor axis.
rmax = 0.4;    % Length of semi-major axis.
Lx = 8;
Ly = 8;
ub = mu*Ly/2/rho/rmax; % shear flow boundary velocity to make Re=1
Grate = 2*ub/Ly;
C = 1/2; %    ds/dx=C fluid-structure grid size constant
Nb = floor(pi*(rmin+rmax)*Nx/Lx/C);
Nb = floor(Nb/4)*4; % membrane grid size

% Time step and final time.
Tfinal = 5.0;
dt = 1e-3; 
NTime = floor(Tfinal./dt)+1;
dt = Tfinal ./ NTime;
tspan = (0:NTime-1)*dt;

angle_true = atan(rmin/rmax*tan(rmin*rmax*Grate*(tspan+dt)/(rmin^2+rmax^2)));

% The initial velocity is zero.
IC_U = @(X,Y)zeros(size(X));
IC_V = @(X,Y)zeros(size(Y));

% The membrane is an ellipse. Align: counter-clockwise from top
IC_ChiX = @(S)Lx/2 - rmin * sin(2*pi*S);
IC_ChiY = @(S)Ly/2 + rmax * cos(2*pi*S);

% The action function.
Action = @( dx, dy, dt, indexT, Lx, Ly, U, V, chiX, chiY, chiX0, chiY0, Fx, Fy)...
    PlotMembrane(Lx, Ly, U, V, chiX,chiY, chiX0, chiY0, 1);

% Do the IB solve.
tic
[X, Y, S, U, V, chiX, chiY, mempos] = ...
        myIBSolver(mu, rho, ub, angle_true, sigma, IC_U, IC_V, IC_ChiX, IC_ChiY,...
        Nx, Ny, Lx, Ly, Nb, NTime, Tfinal, Action);
toc 


%% plot
for indexP = 1:100
    step_l = floor(size(mempos,2)/100);
    i = (indexP-1)*step_l+1;
    chiX = mempos(1:Nb,i);
    chiY = mempos(Nb+1:2*Nb,i);
    chiX0 = mempos(2*Nb+1:3*Nb,i);
    chiY0 = mempos(end-Nb+1:end,i);
    clf;
    axis([0 Lx 0 Ly]);
    xlabel('x'); ylabel('y');

    hold all;

    plot([chiX;chiX(1)],[chiY;chiY(1)],'r');
    plot([chiX0;chiX0(1)],[chiY0;chiY0(1)],'k');


    axis square;
    daspect([1 1 1]);
    drawnow;
    hold off;

    set(gca,'Box','on');
end

angle = -atan((mempos(1,:)-Lx/2)./(mempos(Nb+1,:)-Ly/2));
figure(2)
set(gcf,'position',[300,500,1000,400])
plot(tspan,angle_true,tspan,angle); 
legend('true','numerical');


% Remove PATH reference to avoid clutter
rmpath('./utils');
