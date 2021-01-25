function dydt = fullODE(t,y,Gnew,proj,X,Y,dx,dy,ds,Nb,rho,sigma,L,supp_x,supp_y)
% This function describes a reduced model for the orbit example.
%                du/dt = -Gnew*u+ds/rho*proj*B*f(X)
%                dX/dt = dx*dy*B'*u
% Thus it will be used in MATLAB ODE solvers. e.g. ode45.
% We let y = [u;X], where X=[chiX;chiY].
% The momment matrices M0, M1, M2 are
% 

chiY = y(end-Nb+1:end);
chiX = y(end-2*Nb+1:end-Nb);
u = y(1:end-2*Nb);
B = timeDependentMatB(X,Y,chiX,chiY,dx,dy,supp_x,supp_y);

dXdt = dx*dy*B'*u;

[fx, fy] = elasticity(chiX,chiY,sigma,L,ds,Nb);

dudt = -Gnew*u+ds/rho*proj*B*[fx;fy];
dydt = [dudt;dXdt];