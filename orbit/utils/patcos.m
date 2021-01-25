function [fx,fy] = patcos(x,y)
% Let f = cos(theta_i). f is a function of x_i-1, x_i and x_i+1. This 
% function computes fx1 = df/dx_i+1, fx2 = df/dx_i, fx3 = df/dx_i-1

dx = x([2:end 1])-x;
dy = y([2:end 1])-y;
e = dx.^2+dy.^2;

denom = sqrt(e).*sqrt(e([end 1:end-1]));
num = -dx([end 1:end-1]).*dx-dy([end 1:end-1]).*dy;

fx = ((-dx+dx([end 1:end-1])).*denom-...
        ((-dx+dx([end 1:end-1])).*(-dx).*dx([end 1:end-1])+...
            (-dx).*dy([end 1:end-1]).^2+dx([end 1:end-1]).*dy.^2)./denom.*num)./...
       denom.^2;
   

fy = ((-dy+dy([end 1:end-1])).*denom-...
        ((-dy+dy([end 1:end-1])).*(-dy).*dy([end 1:end-1])+...
            (-dy).*dx([end 1:end-1]).^2+dy([end 1:end-1]).*dx.^2)./denom.*num)./...
       denom.^2;
   