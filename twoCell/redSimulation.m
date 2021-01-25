%% Reduced model simulation in ODE form.
%  two cells interacting with constant background flow


clear;
addpath(genpath('./utils'));

REMAT = true;
RESAMPLE = true;
%% load model parameters
load Model_data/model_params.mat


%% Generate matrix
if REMAT
    fprintf('Computing model matrices...\n');
    [G,proj] = fullMats(Nx,Ny,dx,dy,mu,rho);
    ker = proj*G*proj;
    save Model_data/reduce_model_mat.mat ker G proj
    fprintf('Model matrices saved.\n');
else
    load Model_data/reduce_model_mat.mat 
    fprintf('Loading model matrices finished.\n');
end

%% Sample Moment matrix
supp_x = 4; supp_y = 4;
if RESAMPLE
    fprintf('Sampling moment matrices...\n');
    dxs = .5*dx; dys = .5*dy; %dxs, dys <= 1
    rx = 3*max(rmax1,rmax2); ry = 6*(rmin1+rmin2); 
    [x0,y0,Vs0] = Moment2Dsampling(Xu,Yu,Xv,Yv,dx,dy,Lx,Ly,dxs,dys,proj,rx,ry,supp_x,supp_y);
    [x1,y1,Vs1] = Moment2Dsampling(Xu,Yu,Xv,Yv,dx,dy,Lx,Ly,dxs,dys,ker,rx,ry,supp_x,supp_y);
    save Model_data/sample_MomMat.mat x0 y0 Vs0 x1 y1 Vs1
    fprintf('Sample moment matrices saved.\n');
else
    load Model_data/sample_MomMat.mat 
    fprintf('Loading sample moment matrices finished.\n');
end

%% Simulation
W0 = [U0(:);V0(:)]; W=zeros(2*(Nb1+Nb2),1);

% EULER
Ntime = 1:floor(Tfinal/dt);
tspan = [0 Ntime]*dt;


chi1 = [chiX10;chiY10];
chi1sr = chi1;
chiX1 = chi1(1:end/2); chiX1p = chiX1;
chiY1 = chi1(end/2+1:end); chiY1p = chiY1;
chi2 = [chiX20;chiY20];
chi2sr = chi2;
chiX2 = chi2(1:end/2); chiX2p = chiX2;
chiY2 = chi2(end/2+1:end); chiY2p = chiY2;
chif = [chi1;chi2];
fprintf('Simulation starts...')

for step  = Ntime
tic
    B1 = timeDependentMatB(Xu,Yu,Xv,Yv,chiX1,chiY1,dx,dy,supp_x,supp_y);
    B2 = timeDependentMatB(Xu,Yu,Xv,Yv,chiX2,chiY2,dx,dy,supp_x,supp_y);
    B = [B1 B2];
    
    [fx1,fy1,fx2,fy2] = twoCellForces...
    (sigma1,sigma2,s,ld,a,cd,chiX1,chiY1,chiX10,chiY10,chiX2,chiY2,chiX20,chiY20);
    
    %%% interpolation
    M011 = dataInterp2D(x0,y0,Vs0,chiX1,chiY1,chiX1,chiY1);
    M012 = dataInterp2D(x0,y0,Vs0,chiX1,chiY1,chiX2,chiY2);
    M021 = dataInterp2D(x0,y0,Vs0,chiX2,chiY2,chiX1,chiY1);
    M022 = dataInterp2D(x0,y0,Vs0,chiX2,chiY2,chiX2,chiY2);
    M0 = [M011 M012;M021 M022];
    
    M111 = dataInterp2D(x1,y1,Vs1,chiX1,chiY1,chiX1,chiY1);
    M112 = dataInterp2D(x1,y1,Vs1,chiX1,chiY1,chiX2,chiY2);
    M121 = dataInterp2D(x1,y1,Vs1,chiX2,chiY2,chiX1,chiY1);
    M122 = dataInterp2D(x1,y1,Vs1,chiX2,chiY2,chiX2,chiY2);
    M1 = [M111 M112;M121 M122];
       
    M011p = dataInterp2D(x0,y0,Vs0,chiX1,chiY1,chiX1p,chiY1p);
    M012p = dataInterp2D(x0,y0,Vs0,chiX1,chiY1,chiX2p,chiY2p);
    M021p = dataInterp2D(x0,y0,Vs0,chiX2,chiY2,chiX1p,chiY1p);
    M022p = dataInterp2D(x0,y0,Vs0,chiX2,chiY2,chiX2p,chiY2p);
    M0p = [M011p M012p;M021p M022p];
    M2 = (M0-M0p)/dt;
    
    chiX1p = chiX1; chiX2p = chiX2;
    chiY1p = chiY1; chiY2p = chiY2;
%     %%% Plot
%     Wf = (proj*B*W+W0);
%     U = reshape(Wf(1:Ny*(Nx-1)),Ny,Nx-1);
%     V = reshape(Wf(Ny*(Nx-1)+1:end),Ny-1,Nx);
%     myPlotMembrane(X, Y, U, V,chiX1,chiY1,chiX2,chiY2)
%     %%%
    
    chif = chif+dt*dx*dy*(B'*W0 +M0*W);
    chi1 = chif(1:Nb1*2);
    chi2 = chif(Nb1*2+1:end);
    chi1sr = [chi1sr chi1];
    chi2sr = [chi2sr chi2];
    chiX1 = chi1(1:end/2);
    chiY1 = chi1(end/2+1:end);
    chiX2 = chi2(1:end/2);
    chiY2 = chi2(end/2+1:end);
    
    Mint = -M0\(M1+M2);
    W = W+dt*(Mint*W+1/rho*[ds1*fx1;ds1*fy1;ds2*fx2;ds2*fy2]);
    toc
end

rmpath(genpath('./utils'));

%% plot cells motion
% figure(3)
% axis([3 max(chi1s(end/4+1:end/2,end)) 1.5 4.5])
% daspect([1 1 1])
% for i=1:size(chi1s,2)
%     chi1 = chi1s(:,i);
%     chi2 = chi2s(:,i);
%     chiX1 = chi1(1:end/2);
%     chiY1 = chi1(end/2+1:end);
%     chiX2 = chi2(1:end/2);
%     chiY2 = chi2(end/2+1:end);
%     hold on
% %     plot([chiX1; chiX1(1)],[chiY1; chiY1(1)],'b');
% %     plot([chiX2; chiX2(1)],[chiY2; chiY2(1)],'r');
%     plot(mean(chiX1),mean(chiY1),'b');
%     plot(mean(chiX2),mean(chiY2),'r');
%     drawnow
%     hold off
% end

load Model_data/red_traj
chisr1{cur+1} = chi1sr;
chisr2{cur+1} = chi2sr;
cur = cur + 1;
save Model_data/red_traj chisr1 chisr2 tspan cur

