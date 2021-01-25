%% Reduced projection model simulation.
%  Elliptic membrane with bending force rotates in shear flow with
%  projection model reduced by Krylov subspace method.
% 


clear;
addpath('./utils');

INTERPD = 2; % interpolation dimension
REINTERP = false; % true if regenerate interpolation data
REMAT = false; % true if regenerate model matrices
%% 

%force parameters
sigma = 1e4; % membrane bending coefficient 

rho = 1; % fluid density
mu = 1; % fluid viscosity
ub = 5; % shear flow boundary velocity
Nx = 32; % fluid grid size in x
Ny = 32; % fluid grid size in y
Lx = 2; % fluid range of x
Ly = 2; % fluid range of y
dx = Lx/Nx;
dy = Ly/Ny;
Grate = 2*ub/Ly;

rmin = .1;    % membrane semi-minor axis.
rmax = .15;    % membrane semi-major axis.
C = 1/2; %    ds/dx=C fluid-structure grid size constant
Nb = floor(pi*(rmin+rmax)*Nx/Lx/C);
Nb = floor(Nb/4)*4; % membrane grid size


x = 0:dx:Lx;
y = 0:dy:Ly;

[X,Y] = meshgrid(x,y);

% Initialize fluid velocity and pressure
[U0, u0] = U0_shear(X,Y,ub);
V0 = zeros(size(Y));



% Initialize membrane position
ds = 1/Nb;
S = linspace(0,1-ds,Nb)';
chiX10 = Lx/2+ rmin * sin(2*pi*S);
chiY10 = Ly/2+ rmax * cos(2*pi*S);
chiX01 = chiX10; chiY01 = chiY10;

%% Generate full model matrices
if REMAT
    fprintf('Computing model matrices...\n');
    [G,bvec,proj] = orbMats(Nx,Ny,dx,dy,mu,rho,ub);
    %pGp = proj*G*proj;
    iG = inv(G);
    % Gnew = proj*G;
    % bnew = proj*bvec;
    save Model_data/model_mat.mat G proj iG bvec
    fprintf('Model matrices saved...\n');
else
    fprintf('Loading model matrices...\n');
    load Model_data/model_mat.mat
end
% Gnew = proj*G;
% bnew = proj*bvec;


%% get 2d sample data
if INTERPD == 2
if REINTERP
    fprintf('Computing 2D interpolation data...\n');
    supp_x = 3;
    supp_y = 3;
    blk1 = 1:(Ny-1)*(Nx+1);
    blk2 = (Ny-1)*(Nx+1)+1:size(iG,1);
    dxs = dx/2; dys = dy/2; rx = 5*rmax; ry = 5*rmax; %dxs, dys not too small
    
    iG_blk = iG(blk1,blk1);
    [XiG,YiG,ViG11] = Moment2Dsampling(X,Y,dx,dy,Lx,Ly,dxs,dys,iG_blk,rx,ry,supp_x,supp_y); %sample for iG11
    iG_blk = iG(blk1,blk2);
    [~,~,ViG12] = Moment2Dsampling(X,Y,dx,dy,Lx,Ly,dxs,dys,iG_blk,rx,ry,supp_x,supp_y);%sample for iG12
    iG_blk = iG(blk2,blk1);
    [~,~,ViG21] = Moment2Dsampling(X,Y,dx,dy,Lx,Ly,dxs,dys,iG_blk,rx,ry,supp_x,supp_y); %sample for iG21
    iG_blk = iG(blk2,blk2);
    [~,~,ViG22] = Moment2Dsampling(X,Y,dx,dy,Lx,Ly,dxs,dys,iG_blk,rx,ry,supp_x,supp_y); %sample for iG22
    ViG = {ViG11,ViG12;ViG21,ViG22};

    proj_blk = proj(blk1,blk1);
    [Xpr,Ypr,Vp11] = Moment2Dsampling(X,Y,dx,dy,Lx,Ly,dxs,dys,proj_blk,rx,ry,supp_x,supp_y); %sample for proj11
    proj_blk = proj(blk1,blk2);
    [~,~,Vp12] = Moment2Dsampling(X,Y,dx,dy,Lx,Ly,dxs,dys,proj_blk,rx,ry,supp_x,supp_y); %sample for proj12
    proj_blk = proj(blk2,blk1);
    [~,~,Vp21] = Moment2Dsampling(X,Y,dx,dy,Lx,Ly,dxs,dys,proj_blk,rx,ry,supp_x,supp_y); %sample for proj21
    proj_blk = proj(blk2,blk2);
    [~,~,Vp22] = Moment2Dsampling(X,Y,dx,dy,Lx,Ly,dxs,dys,proj_blk,rx,ry,supp_x,supp_y); %sample for proj22
    Vpr = {Vp11,Vp12;Vp21,Vp22};

    save Model_data/interp_data.mat XiG YiG ViG Xpr Ypr Vpr supp_x supp_y
    fprintf('2D Interpolation data saved...\n');
else
    fprintf('Loading 2D interpolation data...\n');
    load Model_data/interp_data.mat
end
end
%% get 4d sample data
if INTERPD == 4
if REINTERP
    fprintf('Computing 4D interpolation data...\n');
    supp_x = 3;
    supp_y = 3;
    blk1 = 1:(Ny-1)*(Nx+1);
    blk2 = (Ny-1)*(Nx+1)+1:size(iG,1);
    dxs=dx/5; dys=dy/5;

    iG_blk = iG(blk1,blk1);
    [XiG,YiG,ViG11] = Moment4Dsampling(X,Y,dx,dy,dxs,dys,Lx,Ly,iG_blk,2*rmax,2*rmax,supp_x,supp_y); %sample for iG11
    iG_blk = iG(blk1,blk2);
    [~,~,ViG12] = Moment4Dsampling(X,Y,dx,dy,dxs,dys,Lx,Ly,iG_blk,2*rmax,2*rmax,supp_x,supp_y); %sample for iG12
    iG_blk = iG(blk2,blk1);
    [~,~,ViG21] = Moment4Dsampling(X,Y,dx,dy,dxs,dys,Lx,Ly,iG_blk,2*rmax,2*rmax,supp_x,supp_y); %sample for iG21
    iG_blk = iG(blk2,blk2);
    [~,~,ViG22] = Moment4Dsampling(X,Y,dx,dy,dxs,dys,Lx,Ly,iG_blk,2*rmax,2*rmax,supp_x,supp_y); %sample for iG22
    ViG = {ViG11,ViG12;ViG21,ViG22};

    proj_blk = proj(blk1,blk1);
    [Xpr,Ypr,Vp11] = Moment4Dsampling(X,Y,dx,dy,dxs,dys,Lx,Ly,proj_blk,2*rmax,2*rmax,supp_x,supp_y); %sample for proj11
    proj_blk = proj(blk1,blk2);
    [~,~,Vp12] = Moment4Dsampling(X,Y,dx,dy,dxs,dys,Lx,Ly,proj_blk,2*rmax,2*rmax,supp_x,supp_y); %sample for proj12
    proj_blk = proj(blk2,blk1);
    [~,~,Vp21] = Moment4Dsampling(X,Y,dx,dy,dxs,dys,Lx,Ly,proj_blk,2*rmax,2*rmax,supp_x,supp_y); %sample for proj21
    proj_blk = proj(blk2,blk2);
    [~,~,Vp22] = Moment4Dsampling(X,Y,dx,dy,dxs,dys,Lx,Ly,proj_blk,2*rmax,2*rmax,supp_x,supp_y); %sample for proj22
    Vpr = {Vp11,Vp12;Vp21,Vp22};

    save Model_data/interp4d_data.mat XiG YiG ViG Xpr Ypr Vpr supp_x supp_y
    fprintf('4D Interpolation data saved...\n');
else
    fprintf('Loading 4D interpolation data...\n');
    load Model_data/interp4d_data.mat
end
end
%% Initial value
Ub = U0(2:end-1,:); Vb = V0(2:end-1,2:end-1); % more unknown Us bc of Neumann BC (dU/dx=0) at x=0 and x=Lx
W = [Ub(:);Vb(:)];  % original nonzero initial data 
Wf = W;
W1 = zeros(Nb*2,1); % initial data for reduced ODE



%%
Tfinal = 1.4;
dt = 1e-3 ; 
Ntime = 1:floor(Tfinal/dt);
angle = [];

%% simulation using ode solver
odefun = @(t,y) orbredmodel_naive(t,y,iG,proj,X,Y,dx,dy,ds,rho,sigma,...
    chiX01,chiY01,W,4,4); % MaxStep 1e-4
odefun2d = @(t,y) orbredmodel_interp(t,y,X,Y,dx,dy,ds,rho,sigma,...
    chiX01,chiY01,W,supp_x,supp_y,XiG,YiG,ViG,Xpr,Ypr,Vpr);
odefun4d = @(t,y) orbredmodel_interp4D(t,y,X,Y,dx,dy,ds,rho,sigma,...
    chiX01,chiY01,W,supp_x,supp_y,XiG,YiG,ViG,Xpr,Ypr,Vpr);

opts = odeset('MaxStep',2e-5); % for naive only

fprintf('Simulation start...\n');
tic
sol = ode15s(odefun,[0 Tfinal],[W1;chiX10;chiY10],opts);
toc
fprintf('Simulation finished...\n');

tspan = (Ntime-1)*dt;
angle_true = atan(rmin/rmax*tan(rmin*rmax*Grate*tspan/(rmin^2+rmax^2)));

fprintf('Computing rotation angles...\n');
outputs = deval(sol,tspan);

angle = -atan((outputs(end-Nb*2+1,:)-mean(outputs(end-Nb*2+1:end-Nb,:)))./...
    (outputs(end-Nb+1,:)-mean(outputs(end-Nb+1:end,:))));
fprintf('Rotation angle computed...\n');

set(gcf,'position',[300,500,1000,400])
subplot(1,2,1)
chiXend = outputs(end-Nb*2+1:end-Nb,end);
chiYend = outputs(end-Nb+1:end,end);
[chiXendt,chiYendt] = rot_true(chiX01,chiY01,Lx/2,Ly/2,angle_true(end));
plot([chiXend; chiXend(1)],[chiYend; chiYend(1)],'r',[chiXendt; chiXendt(1)],[chiYendt; chiYendt(1)]);
axis([0 Lx 0 Ly]);

subplot(1,2,2)
plot(tspan,angle_true,tspan,angle); 
legend('true','numerical');

rmpath('./utils');
