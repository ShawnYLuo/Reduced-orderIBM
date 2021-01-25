function B = timeDependentMatB_orb(Xu,Yu,Xv,Yv,chiX,chiY,dx,dy,supp_x,supp_y)
chiX = chiX(:)';
chiY = chiY(:)';

xcut = supp_x*dx;
ycut = supp_y*dy;
Bu = phi((Xu(:)-chiX)/xcut,(Yu(:)-chiY)/ycut)/xcut/ycut;
Bv = phi((Xv(:)-chiX)/xcut,(Yv(:)-chiY)/ycut)/xcut/ycut;
% sparsify
B = sparse(blkdiag(Bu,Bv));


