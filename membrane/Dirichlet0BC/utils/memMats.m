function [G,proj,M,N,D2] = memMats(Nx,Ny,dx,dy,mu,rho)
% For the membrane example. The projection model is u'+Gu+Np = Bf, Mu = 0, 
% X' = L^Tu, where G = -mu/rho*[D2,0;0,D2], N = 1/rho*[Dx;Dy], M = [Dx Dy], 
% B and L are generated according to immersed boundary method. With
% projection matrix Sig = N/(M*N)*M and proj = I-Sig, we obtain
%                  u'+proj*G*u = proj*B*f
% where pressure p is eliminated.
%


mysize = Nx*Ny;
% Dx=[0 I 0 ...;-I 0 I 0...;...;0...-I 0] with Nx-by-Nx blocks,
% each block is Ny-by-Ny
vecx = ones((Nx-1)*Ny,1);
Dx = 1/2/dx*(diag(vecx,Ny)-diag(vecx,-Ny));

% Dy
vecy = repmat([ones(Ny-1,1);0],Nx,1);
vecy = vecy(1:end-1);
Dy = 1/2/dy*(diag(vecy,1)-diag(vecy,-1));

% D2 (laplacian)
vec = ones(mysize,1);
Dxx = diag(-2*vec)+diag(vecx,Ny)+diag(vecx,-Ny);
Dyy = diag(-2*vec)+diag(vecy,1)+diag(vecy,-1);
D2 = 1/dx/dx*Dxx+1/dy/dy*Dyy;
%% GMN
G = -mu/rho*[D2 zeros(mysize);zeros(mysize) D2];
M = [Dx Dy];
N = 1/rho*[Dx;Dy];

Sig = N/(M*N)*M;
proj = eye(2*mysize)-Sig;
    