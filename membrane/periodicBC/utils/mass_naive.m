function M1 = mass_naive(t,y,X,Y,dx,dy,Nb,supp_x,supp_y)

chiY = y(end-Nb+1:end);
chiX = y(end-2*Nb+1:end-Nb);
B = timeDependentMatB(X,Y,chiX,chiY,dx,dy,supp_x,supp_y);

M1 = blkdiag(B'*B,eye(2*Nb));