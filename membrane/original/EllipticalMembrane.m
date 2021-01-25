% In this example, we look at an elliptical membrane in a fluid that is
% initially at rest.  The ellipse should oscillate and eventually settle to
% a circle with a prescribed radius.

clear;
% Add PATH reference in order to run solver
addpath('./utils');

% The number of grid points.
N = 2*round(2.^4); 
Nb = 3*N;

% Parameter values.
mu = 1;        % Viscosity.
sigma = 1e3;     % Spring constant.
rho = 1;       % Density.
rmin = 0.2;    % Length of semi-minor axis.
rmax = 0.4;    % Length of semi-major axis.
Lx = 2;
Ly = 2;
L = 1;  % Resting length.

% Time step and final time.
Tfinal = .24;
dt = 1e-4; 
NTime = floor(Tfinal./dt)+1;
dt = Tfinal ./ NTime;

% The initial velocity is zero.
IC_U = @(X,Y)zeros(size(X));
IC_V = @(X,Y)zeros(size(Y));

% The membrane is an ellipse.
IC_ChiX = @(S)Lx/2 + rmax * cos(2*pi*S);
IC_ChiY = @(S)Ly/2 + rmin * sin(2*pi*S);

% The action function.
Action = @( dx, dy, dt, indexT, X, Y, U, V, chiX, chiY, Fx, Fy)...
    PlotMembrane(X, Y, U, V, chiX,chiY,1);

% Do the IB solve.
tic
[X, Y, S, U, V, chiX, chiY, mempos] = ...
        myIBSolver(mu, rho, sigma, L, IC_U, IC_V, IC_ChiX, IC_ChiY,...
        N, N, Lx, Ly, Nb, NTime, Tfinal, Action);
toc 
save ../fftMempos.mat mempos
% Remove PATH reference to avoid clutter
rmpath('./utils');
