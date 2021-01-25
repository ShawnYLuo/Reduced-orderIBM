function dydt = redModel_naive(t,y,iG,proj,Xu,Yu,Xv,Yv,dx,dy,ds,Nb,rho,sigma,rL,supp_x,supp_y)
% This function describes a reduced model for the orbit example.
%             z' = -(M0\M1)*(z+ds/rho*f)
%             X'             = dx*dy*M0*z
% where the matrices
%       M0 = B'*iG*B
%       M1 = B'*proj*B
% These moment matrices are generated naively, i.e. not interpolated and
% thus slower but more accurate.


chiY = y(end-Nb+1:end);
chiX = y(end-2*Nb+1:end-Nb);
z = y(1:end-2*Nb);
B = timeDependentMatB(Xu,Yu,Xv,Yv,chiX,chiY,dx,dy,supp_x,supp_y);
M0 = B'*iG*B;
M1 = B'*proj*B;

dXdt = dx*dy*M0*z;


[fx, fy] = elasticity(chiX,chiY,sigma,rL,ds,Nb);

%Bdot = Bdot_phi(X,Y,chiX,chiY,dXdt(1:end/2),dXdt(end/2+1:end),dx,dy,supp_x,supp_y);


dzdt = (M0\M1)*(-z+ds/rho*[fx;fy]);
dydt = [dzdt;dXdt];