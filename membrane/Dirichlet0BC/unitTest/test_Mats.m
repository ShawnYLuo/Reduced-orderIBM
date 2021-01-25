clear;
Lx = 1;
Ly = 1;
Nx = 6;
Ny = 6;
dx = Lx/Nx;
dy = Ly/Nx;
rho = 1;
mu = 1;

x = linspace(dx,Lx-dx,Nx-1);
y = linspace(dy,Ly-dy,Ny-1);
[X,Y]=meshgrid(x,y);
Z = sin(2*pi/Ly*X).*sin(2*pi/Ly*Y);
mesh(X,Y,Z)

%%
Zx = 2*pi/Lx*cos(2*pi/Lx*X).*sin(2*pi/Ly*Y);
Zy = 2*pi/Lx*sin(2*pi/Ly*X).*cos(2*pi/Ly*Y);
[G,proj,M,N,D2] = memMats(Nx-1,Ny-1,dx,dy,mu,rho);
DZ = rho*N*Z(:);
DZx = DZ(1:end/2);
dx_err = norm(DZx-Zx(:),inf)
DZy = DZ(end/2+1:end);
dy_err = norm(DZy-Zy(:),inf)

subplot(2,2,1)
mesh(X,Y,Zx)
title('true Dx')
subplot(2,2,2)
mesh(X,Y,reshape(DZx,Ny-1,Nx-1))
title('numerical Dx')
subplot(2,2,3)
mesh(X,Y,Zy)
title('true Dy')
subplot(2,2,4)
mesh(X,Y,reshape(DZy,Ny-1,Nx-1))
title('numerical Dy')

