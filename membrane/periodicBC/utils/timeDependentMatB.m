function B = timeDependentMatB(X,Y,chiX,chiY,dx,dy,supp_x,supp_y)
[Ny,Nx] = size(X);
s1=Nx*Ny;
chiX = chiX(:)';
chiY = chiY(:)';
Nb = length(chiX);

xcut = supp_x*dx;
ycut = supp_y*dy;
Bx = phi((X(:)-chiX)/xcut,(Y(:)-chiY)/ycut)/xcut/ycut;
% sparsify
B = sparse([Bx zeros(s1,Nb); zeros(s1,Nb) Bx]);