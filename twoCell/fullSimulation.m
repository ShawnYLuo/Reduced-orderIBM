%% Full model simulation in ODE form.
%  two cells interacting with constant background flow


clear;
addpath(genpath('./utils'));

REMAT = true;
%% load model parameters
load Model_data/model_params.mat

if REMAT
    [G,proj] = fullMats(Nx,Ny,dx,dy,mu,rho);
    Gnew = proj*G;
    save Model_data/fullmodel_mat.mat Gnew proj 
else
    load Model_data/fullmodel_mat.mat 
end

%% Simulation
W0 = [U0(:);V0(:)]; W=zeros(size(W0));

% EULER
Ntime = 1:floor(Tfinal/dt);
tspan = [0 Ntime]*dt;

supp_x = 4; supp_y = 4;
chi1 = [chiX10;chiY10];
chi1sf = chi1;
chiX1 = chi1(1:end/2);
chiY1 = chi1(end/2+1:end);
chi2 = [chiX20;chiY20];
chi2sf = chi2;
chiX2 = chi2(1:end/2);
chiY2 = chi2(end/2+1:end);
fprintf('Simulation starts...')
for step  = Ntime
    Wf = (W+W0);
    tic
    [fx1,fy1,fx2,fy2] = twoCellForces...
    (sigma1,sigma2,s,ld,a,cd,chiX1,chiY1,chiX10,chiY10,chiX2,chiY2,chiX20,chiY20);

%     U = reshape(Wf(1:Ny*(Nx-1)),Ny,Nx-1);
%     V = reshape(Wf(Ny*(Nx-1)+1:end),Ny-1,Nx);
%     myPlotMembrane(X, Y, U, V,chiX1,chiY1,chiX2,chiY2)
    
    B1 = timeDependentMatB(Xu,Yu,Xv,Yv,chiX1,chiY1,dx,dy,supp_x,supp_y);
    B2 = timeDependentMatB(Xu,Yu,Xv,Yv,chiX2,chiY2,dx,dy,supp_x,supp_y);
    chi1= chi1+dt*dx*dy*B1'*Wf;
    chi2= chi2+dt*dx*dy*B2'*Wf;
    chi1sf = [chi1sf chi1];
    chi2sf = [chi2sf chi2];
    W = W+dt*(-Gnew*W+1/rho*proj*(ds1*B1*[fx1;fy1]+ds2*B2*[fx2;fy2]));
    toc
    chiX1 = chi1(1:end/2);
    chiY1 = chi1(end/2+1:end);
    chiX2 = chi2(1:end/2);
    chiY2 = chi2(end/2+1:end);  
end
rmpath(genpath('./utils'));

%% plot cells motion
%figure(3)
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
load Model_data/full_traj
chisf1{cur+1} = chi1sf;
chisf2{cur+1} = chi2sf;
cur = cur + 1;
save Model_data/full_traj chisf1 chisf2 tspan cur
