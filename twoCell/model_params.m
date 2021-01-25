% model parameters that 1) full, 2) bruteforce_red and 3) interpolation_red
% will load.
addpath(genpath('./utils'));
%% force parameters
sigma1 = 1e4; % membrane bending coefficient
sigma2 = 1e4;
a = 1e4; % repulsion force coefficient
s = 1e4; % attraction force coefficient

%% fluid domain
rho = 1; % fluid density
mu = 1; % fluid viscosity
ub = 15; % shear flow boundary velocity
Nx = 64; % fluid grid size in x
Ny = 32; % fluid grid size in y
Lx = 12; % fluid range of x
Ly = 6; % fluid range of y
dx = Lx/Nx;
dy = Ly/Ny;

C = 1/2; %    ds/dx=C fluid-structure grid size constant
rmin1 = .2;    % membrane semi-minor axis.
rmax1 = .3;    % membrane semi-major axis.
Nb1 = floor(pi*(rmin1+rmax1)*Nx/Lx/C);
Nb1 = floor(Nb1/4)*4; % membrane grid size
%Nb1 = 3*Ny;

rmin2 = .2;    % membrane semi-minor axis.
rmax2 = .3;    % membrane semi-major axis.
Nb2 = floor(pi*(rmin2+rmax2)*Nx/Lx/C);
Nb2 = floor(Nb2/4)*4; % membrane grid size
%Nb2 = 3*Ny;

cd = 3*(rmin1+rmin2); % repulsion force separation distance
ld = 3*(rmin1+rmin2); % attraction force separation distance

xf = dx:dx:Lx-dx;
xh = dx/2:dx:Lx-dx/2;
yf = dy:dy:Ly-dy;
yh = dy/2:dy:Ly-dy/2;

[X,Y] = meshgrid(xf,yf);
[Xu,Yu] = meshgrid(xf,yh);
[Xv,Yv] = meshgrid(xh,yf);
[Xp,Yp] = meshgrid(xh,yh);

% Initialize fluid velocity and pressure
%U0 = U0_const(Yu,ub);
U0 = U0_shear(Yu,ub);
V0 = zeros(size(Yv));

% Initialize membrane position
ds1 = 1/Nb1;
S1 = linspace(0,1-ds1,Nb1)';
chiX10 = Lx/3+ rmax1 * sin(2*pi*S1);
chiY10 = (Ly/2-3*rmin1)+ rmin1 * cos(2*pi*S1);

ds2 = 1/Nb2;
S2 = linspace(0,1-ds2,Nb2)';
chiX20 = Lx/3+ rmax2 * sin(2*pi*S2);
chiY20 = (Ly/2+3*rmin2)+ rmin2 * cos(2*pi*S2);

Tfinal = Lx/ub/8;
dt = 5e-4;

save Model_data/model_params.mat sigma1 sigma2 a s cd ld rho mu ub Nx Ny Lx Ly ...
    dx dy rmin1 rmax1 rmin2 rmax2 Nb1 Nb2 X Y Xu Yu Xv Yv Xp Yp U0 V0 ds1...
    chiX10 chiY10 ds2 chiX20 chiY20 Tfinal dt

rmpath(genpath('./utils'));