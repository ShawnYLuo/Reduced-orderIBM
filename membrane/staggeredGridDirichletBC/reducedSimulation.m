%% Reduced projection model simulation.
% TODO: Modify redModel and Bdot accordingly

%  Elliptic membrane with bending force rotates in shear flow with
%  projection model reduced by Krylov subspace method.
% 
%  Hyper parameters: 
%            1) Set REMAT to true to regenerate matrices if model 
%               parameters have been modified.
%            2) Setting USESOLVER to false will run simulation with 
%               forward EULER. If set to true, ode15s will be used for
%               simulation but is not working currently (extremely slow).

clear;
addpath('./utils');

REMAT = false; 
USESOLVER = false;
%% Model parameters
sigma = 1e3; % membrane bending coefficient
rL = 1;  % Resting length.
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
rR = sqrt(rmin*rmax);

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
chiXr = Lx/2+ rR * cos(2*pi*S);
chiYr = Ly/2+ rR * sin(2*pi*S);

%% Generate matrix
if REMAT
    fprintf('Computing model matrices...\n');
    [G,proj,Div,Gradp,Lapu,Lapv] = memMats(Nx,Ny,dx,dy,mu,rho);
    iG = pinv(G);
    save Model_data/reduce_model_mat.mat iG G proj Div Gradp Lapu Lapv
    fprintf('Model matrices saved...\n');
else
    load Model_data/reduce_model_mat.mat 
    fprintf('Loading model matrices finished...\n');
end

%% Simulation
fprintf('Simulation start...\n');
if ~USESOLVER
% Euler
% Simulation parameters
Tfinal = .12;
supp_x = 4; supp_y = 4;
dt = 1e-4; 
Ntime = 1:floor(Tfinal/dt);
tspan = linspace(0,Tfinal,floor(Tfinal/dt)+1);

%z = zeros(4*Nb,1); % 2nd order 
z = zeros(2*Nb,1);
chi = [chiX10;chiY10];

ker = proj*G*proj;

chiX = chi(1:end/2);
chiY = chi(end/2+1:end);
B_prev = timeDependentMatB(Xu,Yu,Xv,Yv,chiX,chiY,dx,dy,supp_x,supp_y);
redMemPos1 = [chiX10(1);chiY10(1)];
tic
for step = Ntime % W = B
    % one step 0.03s
    chiX = chi(1:end/2);
    chiY = chi(end/2+1:end);
    B = timeDependentMatB(Xu,Yu,Xv,Yv,chiX,chiY,dx,dy,supp_x,supp_y);
    dB = (B-B_prev)/dt;
    B_prev = B;

%     % draw fluid and structure    
%     w = proj*B*z;
%     U = reshape(w(1:end/2),Ny,Nx-1);
%     V = reshape(w(end/2+1:end),Ny-1,Nx);
%     myPlotMembrane(X, Y, U, V,chiX,chiY,1)
    
    % each matrix takes about 0.01s to compute
    M0 = B'*proj*B; 
    M1 = B'*ker*B; 
    M2 = B'*proj*dB; 
    
    dchi = dx*dy*M0*z;
    chi = chi+dt*dchi;
    
    [fx, fy] = elasticity(chiX,chiY,sigma,rL,ds,Nb);
    f = [fx;fy];
    z = z+dt*(-M0\(M1+M2)*z+ds/rho*f); % around 0.001s
    redMemPos1 = [redMemPos1 [chiX(1);chiY(1)]];
    %pause
end
toc
fprintf('Simulation finished...\n');
save redMempos1st.mat redMemPos1
load fullMempos.mat
figure(2)
plot(tspan,myMemPos(1,:),tspan,redMemPos1(1,:));
legend('full','reduced1');


rmpath('./utils');

else % simulation using ode solver

% Simulation parameters
Tfinal = .06;
supp_x = 2; supp_y = 2;


odefun = @(t,y) redModel_naive(t,y,iG,proj,Xu,Yu,Xv,Yv,dx,dy,ds,Nb,rho,sigma,rL,supp_x,supp_y); 
%Mass = @(t,y) mass_naive(t,y,X,Y,dx,dy,Nb,supp_x,supp_y);
%opts = odeset('Mass',Mass,'MStateDependence','strong'); 
W1 = zeros(Nb*2,1); % initial data for reduced ODE
tic
sol = ode15s(odefun,[0 Tfinal],[W1;chiX10;chiY10]);
toc
fprintf('Simulation finished...\n');
outputs = deval(sol,tspan);


% plot motion
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

end


