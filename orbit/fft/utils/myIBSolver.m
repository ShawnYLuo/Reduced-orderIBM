function [X, Y, S, U, V, chiX, chiY, mempos] = myIBSolver(mu, rho, ub, angles,...
    sigma, IC_U, IC_V, IC_ChiX, IC_ChiY, Nx, Ny, Lx, Ly, Nb, NTime, Tfinal, ActionFun)
% IBSolver - Solves IB problem
%
% Modification 1) single step method for fluid (myFluidSolver) and 
%                 membrane (update once each time step)
%              2) use different delta function to calculate force 
%                 (myCalculateFoce) and update membrane
%                 (myUpdateMembranePosition)
%              3) track right most point position for comparison
%
% Equation in 2d on a square [0,1]x[0,1]
%    rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) - p_x + line_integral_x
%    rho*v_t = -rho*u*v_x + rho*v*v_y + mu*laplacian(v) - p_y + line_integral_y
%
%    u_x + v_y = 0
%
%    Chi_t = int_\Omega u(x,t) delta(x - Chi(s,t)) dx
%    line_integral_x = int_\Gamma fx(s,t) delta(x - Chi(s,t)) ds
%    line_integral_y = int_\Gamma fy(s,t) delta(x - Chi(s,t)) ds
%
% with Periodic BC, a closed membrane, and initial conditions
%    u(x,y,0) = IC_U(x,y)
%    v(x,y,0) = IC_V(x,y)
%    ChiX(x,y,0) = IC_ChiX(s)
%    ChiY(x,y,0) = IC_ChiY(s)
%
% Parameters:
%   mu                The diffusive constant in navier-stokes.
%   rho               Density of the fluid
%   sigma             Spring constant of the immersed boundary.
%   L                 Resting length of the membrane.
%   IC_U,IC_V,        The initial conditions for the problem.
%   IC_ChiX,IC_ChiY   Make sure they are consistent with the BC. 
%                     Note: IC_ChiX and IC_ChiY are parametrized 
%                     between s=0 to 1.
%   Nx                Number of grid points along the x-axis.
%   Ny                Number of grid points along the y-axis.
%   Lx                Length of domain along the x-axis.
%   Ly                Length of domain along the y-axis.
%   Nb                Number of grid points on membrane
%   NTime             Number of time steps
%   Tfinal            Compute to time
%   ActionFun         Function which allows you do something after 
%                     each timestep (eg. plot the membrane). The function 
%                     needs to have the profile:
%                         ActionAfterTimeStep( dx, dy, dt, indexT, X, Y, u, v, chiX, chiY, Fx, Fy)
%
% Return:
%   X                 A mesh of x grid values.
%   Y                 A mesh of y grid values.
%   S                 A vector of values for s, which parameterises the membrane.
%   U                 The fluid velocity in x direction at the last timestep
%   V                 The fluid velocity in y direction at the last timestep
%   chiX              The membrane x position at the last timestep
%   chiY              The membrane y position at the last timestep
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
% Modified by: Yushuang Luo 7/22/2020

% Nx and Ny need to be even. 
if mod(Nx,2) == 1 && mod(Ny,2) == 1
	error('MATLAB:IBSolver',...
          'The parameters Nx and Ny needs to be even.');
end

% Step Sizes.
% These assume periodicity in both x and s (that is, a closed membrane).
dx = Lx/(Nx);
dy = Ly/(Ny);
dt = Tfinal/NTime;
ds = 1/Nb;  
nTimesVec = 1:NTime;

% Create Mesh.
% These assume periodicity in x, y, and s.
x = (dx:dx:Lx-dx)';
y = (dy:dy:Ly-dy)';
[X,Y] = meshgrid(x,y);
S = linspace(0,1-ds,Nb)';     

% Construct Index Matrices
[IndX,IndY] = meshgrid(0:2*Nx-1,0:2*Ny-1);   

% Construct initial conditions.
U = IC_U(X,Y);
V = IC_V(X,Y);
U0 = repmat(linspace(ub,-ub,Ny+1)',1,Nx+1);
V0 = zeros(size(U0));
chiX = IC_ChiX(S);
chiY = IC_ChiY(S);
chiX0 = chiX;
chiY0 = chiY;

mempos = [];

for indexT=nTimesVec    
    % Step 1: Calculate Force coming from membrane at half time-step
    [Fxh, Fyh, fxh, fyh] = myCalculateForce...
        (chiX, chiY, x, y, ds, dx, dy, Nb, Nx, Ny, sigma, chiX0, chiY0);
   
    % Step 3: Solve for Fluid motion
    [U, V] = myFluidSolve2...
        (U, V, Fxh, Fyh, rho, mu, dx, dy, dt, Nx, Ny, IndX, IndY);

    % Step 4: Update Position of Boundary of membrane again for a half
    % time-step
    [chiX, chiY] = myUpdateMembranePosition...
        (U+U0(2:Ny,2:Nx), V+V0(2:Ny,2:Nx), chiX, chiY, chiX, chiY, x, y, dt, dx, dy, Nb);

    % Call Action Function
    [chiXt, chiYt] = rot_true(chiX0, chiY0, Lx/2, Ly/2, angles(indexT));
    mempos = [mempos [chiX(:);chiY(:); chiXt(:); chiYt(:)]];
     ActionFun( dx, dy, dt, indexT, X, Y, U+U0(2:Ny,2:Nx), V+V0(2:Ny,2:Nx), ...
         chiX, chiY, chiXt, chiYt, Fxh, Fyh);
        
end
