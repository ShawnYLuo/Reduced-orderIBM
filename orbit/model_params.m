% model parameters that 1) fullSimulation and 2) reducedSimulation
% will load.
addpath(genpath('./utils'));
%% force parameters
sigma = 1e4; % membrane bending coefficient

%%
rho = 1; % fluid density
mu = 3; % fluid viscosity
ub = 64/3; % shear flow boundary velocity
Nx = 32; % fluid grid size in x
Ny = 32; % fluid grid size in y
Lx = 8; % fluid range of x
Ly = 8; % fluid range of y
dx = Lx/Nx;
dy = Ly/Ny;
Grate = 2*ub/Ly;

rmin = .2;    % membrane semi-minor axis.
rmax = .3;    % membrane semi-major axis.
C = 1/2; %    ds/dx=C fluid-structure grid size constant
Nb = floor(pi*(rmin+rmax)*Nx/Lx/C);
Nb = floor(Nb/4)*8; % membrane grid size
%Nb = 3*Ny;

xf = dx:dx:Lx-dx;
xh = dx/2:dx:Lx-dx/2;
yf = dy:dy:Ly-dy;
yh = dy/2:dy:Ly-dy/2;

[X,Y] = meshgrid(xf,yf);
[Xu,Yu] = meshgrid(xf,yh);
[Xv,Yv] = meshgrid(xh,yf);
[Xp,Yp] = meshgrid(xh,yh);


% Initialize fluid velocity and pressure
U0 = U0_shear(Yu,ub,Ly);
V0 = zeros(size(Yv));
W0 = [U0(:);V0(:)];


% Initialize membrane position
ds = 1/Nb;
S = linspace(0,1-ds,Nb)';
chiX0 = Lx/2+ rmin * sin(2*pi*S);
chiY0 = Ly/2+ rmax * cos(2*pi*S);

supp_x = 5; supp_y = 5;
one_period = 2*pi*(rmin^2+rmax^2)/rmin/rmax/Grate;
desired_period = 1;
Tfinal = one_period*desired_period;
dt = 2e-4;
save Model_data/model_params.mat sigma rho mu ub Nx Ny Lx Ly ...
    dx dy Grate rmin rmax Nb X Y Xu Yu Xv Yv Xp Yp U0 V0 W0 ds...
    chiX0 chiY0 supp_x supp_y Tfinal dt

rmpath(genpath('./utils'));

