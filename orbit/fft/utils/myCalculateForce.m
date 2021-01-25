function [Fx, Fy, fx, fy] = myCalculateForce...
    (chiX, chiY, x, y, ds, dx, dy, Nb, Nx, Ny, sigma, chiX0, chiY0)
% Modification: use different regularized 2d-delta function phi. The 
% original is the cartesian product of 1d-delta (regularized) functions,
% thus can performa vectorized computation. Not sure how to do it with phi.
% Thus a tiny little bit slower.
%
% This function computes the components of the forcing coming from
% interactions between the fluid and the immersed boundary.
% The components of the force are given by
% F(x,y) = int_Gamma f(s) * delta(x - chiX(s)) * delta(y - chiY(s)) * ds
% where 
% s parameterises the immersed boundary,
% f(s) is the force density on the boundary,
% chiX, chiY are the x and y positions of points on the membrane, 
% delta is a delta function.
%
% The force density is computed using a simple linear spring model with
% resting length L.  This leads to a force density of
%    f = [\chi_s * ( 1 - L / abs(\chi_s) ) ]_s
% where the subscript s refers to a derivative with respect to s, which
% parameterises the membrane.
%
% INPUTS:   chiX        A vector of length Nb that gives the current
%                       x-coordinates of the position of the immersed bounday
%                       at time-step n-1/2.
%           chiY        A vector of length Nb that gives the current
%                       y-coordinates of the position of the immersed bounday
%                       at time-step n-1/2.
%           x,y         Vectors of length N that give the x and y values
%                       occuring in the Eulerian grid.
%           ds          The step size in hte Lagrangian grid.
%           dx          The step size in x-direction for the Eulerian grid.
%           dy          The step size in y-direction for the Eulerian grid.
%           Nb          The total number of grid points in the Lagrangian
%                       grid.
%           Nx          The total number of grid points along x-dimension
%                       of the Eulerian grid.
%           Ny          The total number of grid points along y-dimension
%                       of the Eulerian grid.
%           Lx          The length of the domain along x-direction.
%           Ly          The length of the domain along y-direction.
%           sigma       The spring constant used to compute the force density.
%           L           The resting length of the membrane.
%
% OUTPUTS:  Fx, Fy      The x and y components of the forcing resulting from
%                       the fluid-membrane interactions.  These are N X N
%                       matrices that give the force at each point on the
%                       grid.
%            fx, fy     The x and y components of the force density on the
%                       membrane.  These are Nb X 1 vectors.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

[fx, fy] = bF(sigma, chiX,chiY,chiX0,chiY0);

Fx = zeros(Ny-1,Nx-1);
Fy = zeros(Ny-1,Nx-1);
[X,Y] = meshgrid(x,y);
xcut = 2*dx;
ycut = 2*dy;
for s = 1:Nb
    IntMat = phi((X-chiX(s))/xcut,(Y-chiY(s))/ycut)/xcut/ycut;
    Fx = Fx + fx(s)*ds*IntMat;
    Fy = Fy + fy(s)*ds*IntMat;
end
