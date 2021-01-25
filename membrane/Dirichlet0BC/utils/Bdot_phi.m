function Bdot = Bdot_phi(X,Y,chiX,chiY,dchiX,dchiY,dx,dy,supp_x,supp_y)
% This function calculates dB/dt, where B=[B0 0;0 B0]. Each column of B0 
% is regularized delta function defined on (X,Y) plane. The regularized 
% delta function here is
%               phi(r) = 10/3pi*(1+(2r-3)r^2).
%
[Ny,Nx] = size(X);
s1 = Nx*Ny;
chiX = chiX(:)';
chiY = chiY(:)';
Nb = length(chiX);

Xmat = X(:)-chiX;
Ymat = Y(:)-chiY;

xcut = supp_x*dx;
ycut = supp_y*dy;

% Rijs = sqrt((Xmat/xcut).^2+(Ymat/ycut).^2);


Bdot = -1/xcut/ycut*phi(Xmat/xcut,Ymat/ycut).*(1/xcut/xcut*Xmat*diag(dchiX)+1/ycut/ycut*Ymat*diag(dchiY));

Bdot = [Bdot zeros(s1,Nb); zeros(s1,Nb) Bdot];
Bdot = sparse(Bdot);