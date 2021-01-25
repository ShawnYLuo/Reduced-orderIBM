function [G,proj,M,N,D2] = memMats(Nx,Ny,dx,dy,mu,rho)
% For the membrane example. The projection model is u'+Gu+Np = Bf, Mu = 0, 
% X' = L^Tu, where G = -mu/rho*[D2,0;0,D2], N = 1/rho*[Dx;Dy], M = [Dx Dy], 
% B and L are generated according to immersed boundary method. With
% projection matrix Sig = N/(M*N)*M and proj = I-Sig, we obtain
%                  u'+proj*G*u = proj*B*f
% where pressure p is eliminated.
%


size = Nx*Ny;
% Dx
Dx = zeros(size);
for i = 1:Nx
    l1 = i+1;
    if l1>Nx
        l1 = l1-Nx;
    end
    l2 = i-1;
    if l2<1
        l2 = Nx;
    end
    Dx((i-1)*Ny+1:i*Ny,(l1-1)*Ny+1:l1*Ny) = 1/2/dx*eye(Ny);
    Dx((i-1)*Ny+1:i*Ny,(l2-1)*Ny+1:l2*Ny) = -1/2/dx*eye(Ny);
end

% Dy
Dy = zeros(size);
for i = 1:Nx
    Dy((i-1)*Ny+1:i*Ny,(i-1)*Ny+1:i*Ny) = 1/2/dy*(diag(ones(Ny-1,1),1)+diag(-ones(Ny-1,1),-1)+diag(-1,Ny-1)+diag(1,1-Ny));
end

% D2 (laplacian)
D2 = zeros(size);
for i = 1:Nx
    l1 = i+1;
    if l1>Nx
        l1 = l1-Nx;
    end
    l2 = i-1;
    if l2<1
        l2 = Nx;
    end
    D2((i-1)*Ny+1:i*Ny,(l1-1)*Ny+1:l1*Ny) = 1/dx/dx*eye(Ny);
    D2((i-1)*Ny+1:i*Ny,(l2-1)*Ny+1:l2*Ny) = 1/dx/dx*eye(Ny);
    D2((i-1)*Ny+1:i*Ny,(i-1)*Ny+1:i*Ny) = -2/dx/dx*eye(Ny)+1/dy/dy*(diag(ones(Ny-1,1),1)+diag(ones(Ny-1,1),-1)+diag(1,Ny-1)+diag(1,1-Ny)+diag(-2*ones(Ny,1))); 
end

%% GMN
G = -mu/rho*[D2 zeros(size);zeros(size) D2];
M = [Dx Dy];
N = 1/rho*[Dx;Dy];

Sig = N/(M*N)*M;
proj = eye(2*size)-Sig;
    