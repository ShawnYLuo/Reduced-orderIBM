function dydt = redModel_naive(t,y,G,proj,X,Y,dx,dy,ds,Nb,rho,sigma,rL,supp_x,supp_y)
% This function describes a reduced model for the orbit example.
%             z' = -M2*(M0\z)+ds/rho*M1*f
%             X'             = dx*dy*z
% where the matrices
%       M0 = B'*B
%       M1 = B'*proj*B
%       M2 = B'*proj*G*B-Bdot'*B
% These moment matrices are generated naively, i.e. not interpolated and
% thus slower but more accurate.


chiY = y(end-Nb+1:end);
chiX = y(end-2*Nb+1:end-Nb);
z = y(1:end-2*Nb);
B = timeDependentMatB(X,Y,chiX,chiY,dx,dy,supp_x,supp_y);
M0 = B'*B;
M1 = B'*proj*B;

dXdt = dx*dy*z;


[fx, fy] = elasticity(chiX,chiY,sigma,rL,ds,Nb);

Bdot = Bdot_phi(X,Y,chiX,chiY,dXdt(1:end/2),dXdt(end/2+1:end),dx,dy,supp_x,supp_y);
M2 = B'*proj*G*B-Bdot'*B;

dzdt = -M2*(M0\z)+ds/rho*M1*[fx;fy];
dydt = [dzdt;dXdt];