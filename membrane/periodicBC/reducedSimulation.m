%% Reduced projection model simulation.
%  Elliptic membrane with bending force rotates in shear flow with
%  projection model reduced by Krylov subspace method.
% 


clear;
addpath('./utils');

REMAT = true; % true if regenerate model matrices
%% force parameters
sigma = 1e2; % membrane bending coefficient
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


x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;

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
    fprintf('Computing model matrices...\n');
    [G,proj,~,~,~] = memMats(Nx,Ny,dx,dy,mu,rho);
    iG = pinv(G);
    save Model_data/reduce_model_mat.mat iG proj G
    fprintf('Model matrices saved...\n');
else
    fprintf('Loading model matrices...\n');
    load Model_data/reduce_model_mat.mat 
end



%% simulation using ode solver
% Simulation parameters
Tfinal = .01;
supp_x = 2; supp_y = 2;


odefun = @(t,y) redModel_naive(t,y,G,proj,X,Y,dx,dy,ds,Nb,rho,sigma,...
    rL,supp_x,supp_y); 
%Mass = @(t,y) mass_naive(t,y,X,Y,dx,dy,Nb,supp_x,supp_y);
%opts = odeset('Mass',Mass,'MStateDependence','strong'); 
W1 = zeros(Nb*2,1); % initial data for reduced ODE
fprintf('Simulation start...\n');
tic
sol = ode15s(odefun,[0 Tfinal],[W1;chiX10;chiY10]);
toc
fprintf('Simulation finished...\n');
outputs = deval(sol,tspan);
%% plot simulation

dt = 1e-4; 
Ntime = 1:floor(Tfinal/dt);
tspan = [0, Ntime]*dt;
redMemPos = [];
for step = 1:size(outputs,2)
    chiX = outputs(end-2*Nb+1:end-Nb,step);
    chiY = outputs(end-Nb+1:end,step);
    
    redMemPos = [redMemPos [chiX(1);chiY(1)]];
    
    plot([chiX; chiX(1)],[chiY; chiY(1)],'r');
    axis([0 Lx 0 Ly]);
    daspect([1 1 1])
end
load fullMempos.mat
figure(2)
plot(tspan,myMemPos(1,1:length(tspan)),tspan,redMemPos(1,:));
legend('full','reduced');


rmpath('./utils');
return


%% Euler
% Simulation parameters
Tfinal = .04;
supp_x = 2; supp_y = 2;
dt = 1e-4 ; 
Ntime = 1:floor(Tfinal/dt);
tspan = [0, Ntime]*dt;
redMemPos = [];


Kernel = proj*G*proj;

z = zeros(2*Nb,1);
chi = [chiX10;chiY10];
for step = Ntime
    chiX = chi(1:end/2);
    chiY = chi(end/2+1:end);
    redMemPos = [redMemPos [chiX(1);chiY(1)]];
    
    B = timeDependentMatB(X,Y,chiX,chiY,dx,dy,supp_x,supp_y);
    M0 = B'*proj*B;
    M1 = B'*Kernel*B;
    [fx, fy] = elasticity(chiX,chiY,sigma,rL,ds,Nb);
    dchi = dx*dy*z;
    chi = chi+dt*dchi;
    Bdot = Bdot_phi(X,Y,chiX,chiY,dchi(1:end/2),dchi(end/2+1:end),dx,dy,supp_x,supp_y);
    M2 = Bdot'*proj*B;
    
    
    z = z+dt*((M2-M1)*pinv(M0)*z+ds/rho*M0*[fx;fy]);
    w = proj*B*pinv(M0)*z;
    U = reshape(w(1:Ny*Nx),Ny,Nx);
    V = reshape(w(1+Ny*Nx:2*Ny*Nx),Ny,Nx);
    clf
    axis([0 Lx 0 Ly]);
    daspect([1 1 1]);
    hold on
    quiver(X,Y,U,V);
    plot([chiX; chiX(1)],[chiY; chiY(1)],'r');
    drawnow;
    hold off
    %pause
end

load fullMempos.mat
figure(2)
plot(tspan,myMemPos(1,1:length(tspan)),tspan,redMemPos(1,:));
legend('full','reduced');


rmpath('./utils');