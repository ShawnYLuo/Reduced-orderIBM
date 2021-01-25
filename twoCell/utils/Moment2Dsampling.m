function [x,y,Vs] = Moment2Dsampling(Xu,Yu,Xv,Yv,dx,dy,Lx,Ly,dxs,dys,M,rx,ry,supp_x,supp_y)
% This function samples the 2D function that approximately generates the
% entries of a block of the moment matrix M=[M11 M12;M21 M22].
% Inputs:  X,Y:      Fluid domain meshgrid
%          dx,dy:    Fluid grid size
%          Lx,Ly:    Fluid domain size
%          dxs,dys:  Sample points grid size
%          Mij:      a block of M
%          rx,ry:    range of sample points
% Outputs: Xs,Ys:    sample points meshgrid
%          Vs:       sample points values
%

% construct sample points mesh
su = length(Xu(:));
Nxs = floor(rx/dxs)*2+1; Nys = floor(ry/dys)*2+1;

x = linspace(-rx,rx,Nxs);
y = linspace(-ry,ry,Nys);
chiXl = Lx/2*ones(size(y));
chiYl = Ly/2+y;

chiXr = Lx/2-x;
chiYr = Ly/2*ones(size(x));

xcut = supp_x*dx;
ycut = supp_y*dy;
Blu = phi((Xu(:)-chiXl(:)')/xcut,(Yu(:)-chiYl(:)')/ycut)/xcut/ycut;
Blv = phi((Xv(:)-chiXl(:)')/xcut,(Yv(:)-chiYl(:)')/ycut)/xcut/ycut;

Bru = phi((Xu(:)-chiXr(:)')/xcut,(Yu(:)-chiYr(:)')/ycut)/xcut/ycut;
Brv = phi((Xv(:)-chiXr(:)')/xcut,(Yv(:)-chiYr(:)')/ycut)/xcut/ycut;

Vs = cell(2,2);
M11 = M(1:su,1:su); Vs{1,1} = Blu'*M11*Bru;
M12 = M(1:su,su+1:end); Vs{1,2} = Blu'*M12*Brv;
M21 = M(su+1:end,1:su); Vs{2,1} = Blv'*M21*Bru;
M22 = M(su+1:end,su+1:end); Vs{2,2} = Blv'*M22*Brv;
