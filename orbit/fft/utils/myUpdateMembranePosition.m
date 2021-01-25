 function [chiX, chiY] = myUpdateMembranePosition...
     (U, V, chiX, chiY, chiXConv, chiYConv, x, y, dt, dx, dy, Nb)
% 
% This function updates the position of the immersed boundary by evolving
% the equation 
%     Chi_t = int_\Omega u(x,t) delta(x - Chi(s,t)) dx.
%
% INPUTS:  U,V                 The x and y components of the fluid velocity at the
%                              previous time step, in matrix format.
%          chiX, chiY          The x and y component of the membrane position at the
%                              previous time step n-1.  These are vectors.
%          chiXConv, chiYConv  The x and y component of the membrane position
%                              used in the integral (summation). These are vectors.                     
%          x, y                The x and y positions in the grid.  These are vectors
%                              of length N.
%          dt                  The time step, a scalar.
%          dx                  The grid spacing, a scalar, in x-direction.
%          dy                  The grid spacing, a scalar, in y-direction.
%          Nb                  The total number of points used to describe the
%                              membrane position.
%          Nx                  The total number of grid points along x-direction
%                              of the domain.
%          Ny                  The total number of grid points along y-direction
%                              of the domain.
%          Lx                  The length of the domain in x-direction.
%          Ly                  The length of the domain in y-direction.
%
% OUTPUTS: chiX, chiY The updated x and y components of the membrane
%                     position.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

[X,Y] = meshgrid(x,y);
xcut = 2*dx;
ycut = 2*dy;
for s = 1:Nb
    IntMat = phi((X-chiXConv(s))/xcut,(Y-chiYConv(s))/ycut)/xcut/ycut;
    chiX(s) = chiX(s) + dt*dx*dy*sum(U.*IntMat,'all');
    chiY(s) = chiY(s) + dt*dx*dy*sum(V.*IntMat,'all');
end