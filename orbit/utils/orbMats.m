function [G,proj] = orbMats(Nx,Ny,dx,dy,mu,rho)
% For the orbit example, the projection model is u'+Gu+Np = Bf, Mu = 0, 
% X' = L^T(u+u0), where G = -mu/rho*[D2,0;0,D2], N = 1/rho*[Dx;Dy], M = [Dx Dy],
% We use 0 Dirichlet boundary conditions for v. For boundary of u, we use 
% 0 Neumann on left and right, and 0 Dirichlet on upper and bottom. 
% B and L are generated according to immersed boundary method. 
% With projection matrix Sig = N/(M*N)*M and proj = I-Sig, we obtain
%                  u'+proj*G*u = proj*B*f
% where pressure p is eliminated. Note that the background flow is not 
% included.
%


% Dx=[-I I 0...;0 -I I 0...;...;0...0 -I I] with (Nx-1)-by-Nx blocks,
% each block is Ny-by-Ny
Dx = 1/dx*(diag(-ones(Nx*Ny,1))+diag(ones(Nx*Ny-Ny,1),Ny));
Dx = Dx(1:Nx*Ny-Ny,:);
% Dy
Bl = -diag(ones(Ny,1))+diag(ones(Ny-1,1),1);
Bl = Bl(1:Ny-1,:);
Dy = Bl;
for i=1:Nx-1
    Dy = blkdiag(Dy,Bl);
end
Dy = 1/dy*Dy;
% D2 (laplacian)
vecu = repmat([ones(Ny-1,1);0],Nx-1,1);
Dxx = diag(-2*ones(Nx*Ny-Ny,1))+diag(ones((Nx-2)*Ny,1),Ny)+diag(ones((Nx-2)*Ny,1),-Ny);
Dxx(1:Ny,1:Ny) = Dxx(1:Ny,1:Ny) + eye(Ny);
Dxx(end-Ny+1:end,end-Ny+1:end) = Dxx(end-Ny+1:end,end-Ny+1:end) + eye(Ny);
Dyy = diag(repmat([-3;-2*ones(Ny-2,1);-3],Nx-1,1))+diag(vecu(1:end-1),1)+diag(vecu(1:end-1),-1);
D2u = 1/dx/dx*Dxx+1/dy/dy*Dyy;

vecv = repmat([ones(Ny-2,1);0],Nx,1);
Dxx = diag([-3*ones(Ny-1,1);-2*ones((Nx-2)*(Ny-1),1);-3*ones(Ny-1,1)])+...
    diag(ones((Nx-1)*(Ny-1),1),Ny-1)+diag(ones((Nx-1)*(Ny-1),1),-Ny+1);
Dyy = diag(-2*ones(Nx*Ny-Nx,1))+diag(vecv(1:end-1),1)+diag(vecv(1:end-1),-1);
D2v = 1/dx/dx*Dxx+1/dy/dy*Dyy;
%% GMN
G = -mu/rho*blkdiag(D2u,D2v);
M = [Dx' Dy'];
M(1:Ny,1:(Nx-1)*Ny)=zeros(Ny,(Nx-1)*Ny);
M(end-Ny+1:end,1:(Nx-1)*Ny)=zeros(Ny,(Nx-1)*Ny);
N = 1/rho*[Dx;Dy];

Sig = N*pinv(M*N)*M;
proj = eye(size(Sig))-Sig;