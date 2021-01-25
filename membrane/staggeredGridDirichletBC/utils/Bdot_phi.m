function Bdot = Bdot_phi(Xu,Yu,Xv,Yv,chiX,chiY,dchiX,dchiY,dx,dy,supp_x,supp_y)
% This function calculates dB/dt, where B=[Bu 0;0 Bv]. Each column of Bu 
% (or Bv) is regularized delta function defined on (X,Y) plane. 
% The regularized delta function here is
%               phi(r) = 10/3pi*(1+(2r-3)r^2).
%


xcut = supp_x*dx;
ycut = supp_y*dy;
% Bdot_u
chiX = chiX(:)';
chiY = chiY(:)';

Xmat_u = Xu(:)-chiX;
Ymat_v = Yu(:)-chiY;

% Rijs = sqrt((Xmat/xcut).^2+(Ymat/ycut).^2);

Bdot_u = -1/xcut/ycut*phi(Xmat_u/xcut,Ymat_v/ycut).*(1/xcut/xcut*Xmat_u*diag(dchiX)+1/ycut/ycut*Ymat_v*diag(dchiY));


%Bdot_v
chiX = chiX(:)';
chiY = chiY(:)';

Xmat_v = Xv(:)-chiX;
Ymat_v = Yv(:)-chiY;



Bdot_v = -1/xcut/ycut*phi(Xmat_v/xcut,Ymat_v/ycut).*(1/xcut/xcut*Xmat_v*diag(dchiX)+1/ycut/ycut*Ymat_v*diag(dchiY));

Bdot = sparse(blkdiag(Bdot_u,Bdot_v)); 