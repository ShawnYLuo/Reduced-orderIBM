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
rL = 0;  % Resting length.
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


x = dx:dx:Lx-dx;
y = dy:dy:Ly-dy;

[X,Y] = meshgrid(x,y);

% Initialize fluid velocity
U0 = zeros(size(X));
V0 = zeros(size(Y));
W0 = [U0(:);V0(:)];


% Initialize membrane position
ds = 1/Nb;
S = linspace(0,1-ds,Nb)';
chiX10 = Lx/2+ rmax * cos(2*pi*S);
chiY10 = Ly/2+ rmin * sin(2*pi*S);

%% Generate matrix
if REMAT
    [G,proj,~,~,~] = memMats(Nx-1,Ny-1,dx,dy,mu,rho);
    Gnew = proj*G;
    save Model_data/fullmodel_mat.mat Gnew proj 
else
    load Model_data/fullmodel_mat.mat 
end

%% Simulation 
if ~USESOLVER
% EULER
Tfinal = .06;
supp_x = 2; supp_y = 2;
dt = 1e-4;
chi = [chiX10;chiY10];
myMemPos = [chiX10(1);chiY10(1)];
tspan = linspace(0,Tfinal,floor(Tfinal/dt)+1);
for step  = 1:floor(Tfinal/dt)
    U = reshape(W0(1:end/2),Ny-1,Nx-1);
    V = reshape(W0(end/2+1:end),Ny-1,Nx-1);
    myPlotMembrane(X, Y, U, V,chiX10,chiY10,1)
    
    [fx,fy] = elasticity(chiX10,chiY10,sigma,rL,ds,Nb);
    B = timeDependentMatB(X,Y,chiX10,chiY10,dx,dy,supp_x,supp_y);
    W0 = W0+dt*(-Gnew*W0+ds/rho*proj*B*[fx;fy]);
    chi= chi+dt*dx*dy*B'*W0;
    chiX10 = chi(1:end/2);
    chiY10 = chi(end/2+1:end);
    myMemPos = [myMemPos [chiX10(1);chiY10(1)]];
end

load ../fftMempos.mat
figure(2)
plot(tspan,mempos(1,1:length(tspan)),tspan,myMemPos(1,:));
legend('fft','fullODE');
save fullMempos.mat myMemPos 


else
% Simulation by ODE solver
Tfinal = .06;

supp_x = 2; supp_y=2;
odefun = @(t,y) fullODE(t,y,Gnew,proj,X,Y,dx,dy,ds,Nb,rho,sigma,rL,supp_x,supp_y);
tic
sol = ode45(odefun,[0 Tfinal],[W0;chiX10;chiY10]);
toc
% plot simulation
dt = 1e-4 ; 
tspan = linspace(0,Tfinal,floor(Tfinal/dt));
outputs = deval(sol,tspan);

myMemPos = [];
for step = 1:size(outputs,2)
    U = reshape(outputs(1:(Ny-1)*(Nx-1),step),Ny-1,Nx-1);
    V = reshape(outputs(1+(Ny-1)*(Nx-1):end-2*Nb,step),Ny-1,Nx-1);
    chiX = outputs(end-2*Nb+1:end-Nb,step);
    chiY = outputs(end-Nb+1:end,step);
    
    myMemPos = [myMemPos [chiX(1);chiY(1)]];
    myPlotMembrane(X, Y, U, V,chiX,chiY,1)
end

load ../fftMempos.mat
figure(2)
plot(tspan,mempos(1,1:length(tspan)),tspan,myMemPos(1,:));
legend('fft','fullODE');
save fullMempos.mat myMemPos
end
rmpath('./utils');
