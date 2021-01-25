%% Full projection model simulation.
%  Elastic membrane in stationary flow with elasticity with projection
%  model instead of ACM.
%  


clear;
addpath('./utils');
REMAT = false; % if regenerate projection and G
USESOLVER = false;
%% force parameters
sigma = 1e3; % membrane bending coefficient
rL = 1;  % Resting length.
%%
rho = 1; % fluid density
mu = 1; % fluid viscosity
Nx = 32; % fluid grid size in x
Ny = 32; % fluid grid size in y
Lx = 2; % fluid range of x
Ly = 2; % fluid range of y
dx = Lx/Nx;
dy = Ly/Ny;


rmin = .2;    % membrane semi-minor axis.
rmax = .4;    % membrane semi-major axis.

Nb = 3*Nx;


xf = dx:dx:Lx-dx;
xh = dx/2:dx:Lx-dx/2;
yf = dy:dy:Ly-dy;
yh = dy/2:dy:Ly-dy/2;

[X,Y] = meshgrid(xf,yf);
[Xu,Yu] = meshgrid(xf,yh);
[Xv,Yv] = meshgrid(xh,yf);
[Xp,Yp] = meshgrid(xh,yh);

% Initialize fluid velocity
U0 = zeros(size(Xu));
V0 = zeros(size(Yv));
W0 = [U0(:);V0(:)];


% Initialize membrane position
ds = 1/Nb;
S = linspace(0,1-ds,Nb)';
chiX10 = Lx/2+ rmax * cos(2*pi*S);
chiY10 = Ly/2+ rmin * sin(2*pi*S);

%% Generate matrix
if REMAT
    [G,proj,~,~,~,~] = memMats(Nx,Ny,dx,dy,mu,rho);
    Gnew = proj*G;
    save Model_data/fullmodel_mat.mat Gnew proj 
else
    load Model_data/fullmodel_mat.mat 
end

%% Simulation 
if ~USESOLVER
% EULER
Tfinal = .12;
supp_x = 4; supp_y = 4;
dt = 1e-4;
chi = [chiX10;chiY10];
chiX = chi(1:end/2);
chiY = chi(end/2+1:end);
myMemPos = [chiX(1);chiY(1)];
tspan = linspace(0,Tfinal,1201);
tic
for step  = 1:1200
    % one step 0.3s
%     U = reshape(W0(1:end/2),Ny,Nx-1);
%     V = reshape(W0(end/2+1:end),Ny-1,Nx);
%     myPlotMembrane(X, Y, U, V,chiX,chiY,1)

    [fx,fy] = elasticity(chiX,chiY,sigma,rL,ds,Nb);
    B = timeDependentMatB(Xu,Yu,Xv,Yv,chiX,chiY,dx,dy,supp_x,supp_y);
    W0 = W0+dt*(-Gnew*W0+ds/rho*proj*B*[fx;fy]);
    chi= chi+dt*dx*dy*B'*W0;
    chiX = chi(1:end/2);
    chiY = chi(end/2+1:end);    
    myMemPos = [myMemPos [chiX(1);chiY(1)]];
end
toc
save fullMempos.mat myMemPos 
save full_eq.mat chi  W0
load ../fftMempos.mat
figure(2)
plot(tspan,mempos(1,:),tspan,myMemPos(1,:));
legend('fft','fullODE');



else
% Simulation by ODE solver
Tfinal = .24;

supp_x = 2; supp_y=2;
odefun = @(t,y) fullODE(t,y,Gnew,proj,Xu,Yu,Xv,Yv,dx,dy,ds,Nb,rho,sigma,rL,supp_x,supp_y);
tic
sol = ode45(odefun,[0 Tfinal],[W0;chiX10;chiY10]);
toc
% plot simulation
dt = 1e-4 ; 
tspan = linspace(0,Tfinal,floor(Tfinal/dt));
outputs = deval(sol,tspan);

myMemPos = [];
for step = 1:size(outputs,2)
    U = reshape(outputs(1:Nx*Ny-Ny,step),Ny,Nx-1);
    V = reshape(outputs(1+Ny*(Nx-1):end-2*Nb,step),Ny-1,Nx);
    chiX = outputs(end-2*Nb+1:end-Nb,step);
    chiY = outputs(end-Nb+1:end,step);
    
    myMemPos = [myMemPos [chiX(1);chiY(1)]];
    myPlotMembrane(X, Y, U, V,chiX,chiY,1);
end

load ../fftMempos.mat
figure(2)
plot(tspan,mempos(1,1:length(tspan)),tspan,myMemPos(1,:));
legend('fft','fullODE');
save fullMempos.mat myMemPos
end
rmpath('./utils');
