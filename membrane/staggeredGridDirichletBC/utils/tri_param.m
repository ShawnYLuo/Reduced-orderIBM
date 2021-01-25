function [a,b] = tri_param(x,y)
% The membrane position is X_i = (x_i,y_i), i =1,2,...,N. This
% function computes the sine and cosine values of the angle between edges.
% b_i = cos(theta_i) and a_i = sin(theta_i), where
% theta_i=angle(X_i-1,X_i,X_i+1).

%compute edge length
dx = x([2:end 1])-x;
dy = y([2:end 1])-y;
e = dx.^2+dy.^2;


dx2 = x([2:end 1])-x([end 1:end-1]);
dy2 = y([2:end 1])-y([end 1:end-1]);
c = dx2.^2+dy2.^2;

k = norm(c);
e = e/k;
c = c/k;

b = (e+e([end 1:end-1])-c)./(2*sqrt(e.*e([end 1:end-1])));
a = sqrt(1-b.^2);