%% Full model simulation in ODE form.
%  Elliptic membrane with bending force rotates in shear flow. 


clear;
addpath('./utils');

USESOLVER = false;
REMAT = true;
%% load model parameters
load Model_data/model_params.mat

%% Generate matrix
if REMAT
    [G,proj] = orbMats(Nx,Ny,dx,dy,mu,rho);
    Gnew = proj*G;
    save Model_data/fullmodel_mat.mat Gnew proj 
else
    load Model_data/fullmodel_mat.mat 
end

%% Simulation
W=zeros(size(W0));

if ~USESOLVER
% EULER
Ntime = 1:floor(Tfinal/dt);
tspan = [0 Ntime]*dt;
angle_true = atan(rmin/rmax*tan(rmin*rmax*Grate*tspan/(rmin^2+rmax^2)));

angle = 0;

chi = [chiX0;chiY0];
memlocs = chi;
chiX = chi(1:end/2);
chiY = chi(end/2+1:end);
fprintf('Simulation starts...')

for step  = Ntime
    %tic
     Wf = (W+W0);
%     U = reshape(Wf(1:Ny*(Nx-1)),Ny,Nx-1);
%     V = reshape(Wf(Ny*(Nx-1)+1:end),Ny-1,Nx);
%     [chiXt, chiYt] = rot_true(chiX0,chiY0,Lx/2,Ly/2,angle_true(step));
%     myPlotMembrane(X, Y, U, V,chiX,chiY,chiXt,chiYt,1)

    [fx, fy] = bF(sigma,chiX,chiY,chiX0,chiY0);
    B = timeDependentMatB_orb(Xu,Yu,Xv,Yv,chiX,chiY,dx,dy,supp_x,supp_y);
    chi= chi+dt*dx*dy*B'*Wf;
    W = W+dt*(-Gnew*W+ds/rho*proj*B*[fx;fy]);
    chiX = chi(1:end/2);
    chiY = chi(end/2+1:end);   
    angle = [angle atan((chiX(Nb/2+1)-chiX(1))/(chiY(1)-chiY(Nb/2+1)))];
    memlocs = [memlocs chi];
    %toc
end
%%
%angle = atan((memlocs(Nb/2+1,:)-memlocs(1,:))./(memlocs(Nb+1,:)-memlocs(Nb/2*3+1)));
%angle = atan((Lx/2-memlocs(1,:))./(memlocs(Nb+1,:)-Ly/2));
figure(2)
plot(tspan,angle_true,tspan,angle); 
legend('true','full');
fullangle = angle;
fullmemlocs = memlocs;
save fullangle.mat fullangle fullmemlocs
else
% simulation using ode45
Tfinal = 1.4;
dt = 1e-4 ; % 
Ntime = 1:floor(Tfinal/dt);
angle = [];

supp_x = 3; supp_y=3;
odefun = @(t,y) orbfullmodel(t,y,Gnew,proj,X,Y,dx,dy,ds,rho,sigma,chiX01,chiY01,W0,supp_x,supp_y);
tic
sol = ode45(odefun,[0 Tfinal],[W1;chiX10;chiY10]);
toc
tspan = (Ntime-1)*dt;
angle_true = atan(rmin/rmax*tan(rmin*rmax*Grate*tspan/(rmin^2+rmax^2)));


outputs = deval(sol,tspan);

angle = -atan((outputs(end-Nb*2+1,:)-mean(outputs(end-Nb*2+1:end-Nb,:)))./...
    (outputs(end-Nb+1,:)-mean(outputs(end-Nb+1:end,:))));

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

end
rmpath('./utils');