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

REMAT = true; 
RESAMPLE = true;

%% load model parameters
load Model_data/model_params.mat


%% Generate matrix
if REMAT
    fprintf('Computing model matrices...\n');
    [G,proj] = orbMats(Nx,Ny,dx,dy,mu,rho);
    ker = proj*G*proj;
    save Model_data/reduce_model_mat.mat ker G proj
    fprintf('Model matrices saved.\n');
else
    load Model_data/reduce_model_mat.mat 
    fprintf('Loading model matrices finished.\n');
end

%% Sample Moment matrix
if RESAMPLE
    fprintf('Sampling moment matrices...\n');
    dxs = dx; dys = dy; %dxs, dys <= 1
    rx = 3*rmax; ry = 3*rmax; 
    [x0,y0,Vs0] = Moment2Dsampling(Xu,Yu,Xv,Yv,dx,dy,Lx,Ly,dxs,dys,proj,rx,ry,supp_x,supp_y);
    [x1,y1,Vs1] = Moment2Dsampling(Xu,Yu,Xv,Yv,dx,dy,Lx,Ly,dxs,dys,ker,rx,ry,supp_x,supp_y);
    save Model_data/sample_MomMat.mat x0 y0 Vs0 x1 y1 Vs1
    fprintf('Sample moment matrices saved.\n');
else
    load Model_data/sample_MomMat.mat 
    fprintf('Loading sample moment matrices finished.\n');
end

%% Simulation
fprintf('Simulation start...\n');
% Euler
% Simulation parameters 
Ntime = 1:floor(Tfinal/dt);
tspan = [0 Ntime]*dt;
angle_true = atan(rmin/rmax*tan(rmin*rmax*Grate*tspan/(rmin^2+rmax^2)));

angle = 0;

chii = [chiX0;chiY0];
memlocs = chii;
chiXi = chii(1:end/2);
chiYi = chii(end/2+1:end);

zi = zeros(2*Nb,1);

chiX_prev = chii(1:end/2);
chiY_prev = chii(end/2+1:end);
for step = Ntime 
    %tic
    B_prev = timeDependentMatB_orb(Xu,Yu,Xv,Yv,chiXi,chiYi,dx,dy,supp_x,supp_y);
%     Wf = proj*B_prev*zi+W0;
%     U = reshape(Wf(1:Ny*(Nx-1)),Ny,Nx-1);
%     V = reshape(Wf(Ny*(Nx-1)+1:end),Ny-1,Nx);
%     myPlotMembrane(X, Y, U, V,chiXi,chiYi,1)
	    
    M0i = dataInterp2D(x0,y0,Vs0,chiXi,chiYi,chiXi,chiYi);
    M1i = dataInterp2D(x1,y1,Vs1,chiXi,chiYi,chiXi,chiYi);
    %M0d = dataInterp2D(x0,y0,Vs0,chiXi,chiYi,chiX_prev,chiY_prev);
    %M2i = (M0i-M0d)/dt;
    
    chiX_prev = chiXi; chiY_prev = chiYi;
	
	[fxi, fyi] = bF(sigma,chiXi,chiYi,chiX0,chiY0);
    fi = [fxi;fyi];
    
	dchii = dx*dy*(M0i*zi+B_prev'*W0);
    chii = chii+dt*dchii;
	chiXi = chii(1:end/2);
    chiYi = chii(end/2+1:end);
    
    %Mint = -M0i\(M1i+M2i);
    Mint = -M0i\M1i;
    zi = zi+dt*(Mint*zi+ds/rho*fi); % around 0.001s
	
	angle = [angle atan((chiXi(Nb/2+1)-chiXi(1))/(chiYi(1)-chiYi(Nb/2+1)))];
    memlocs = [memlocs chii];
    %toc
end
fprintf('Simulation finished...\n');
angle_true = atan(rmin/rmax*tan(rmin*rmax*Grate*tspan/(rmin^2+rmax^2)));
%%
figure(2)
load fullangle
plot(tspan,angle_true,'k',tspan,fullangle,'b--',tspan,angle,':r','LineWidth',3);
axis([0 2.6 -2 2])
legend('Analtyical','FOM','ROM','Orientation','horizontal');



rmpath('./utils');




