function [U, V] = myFluidSolve(U, V, Fx, Fy, rho, mu, dx, dy, dt, Nx, Ny, Lx, Ly, IndX, IndY)
% Modification: no advection thus no need to do two-step
%
% This function solves the incompressible Navier-Stokes (NS) equations
%      rho*u_t = -(rho*u*u_x + rho*v*u_y) + mu*laplacian(u) - p_x + Fx
%      rho*v_t = -(rho*u*v_x + rho*v*v_y) + mu*laplacian(v) - p_y + Fy
%      u_x + v_y = 0.
%
% INPUTS:  U, V           The x and y components of the fluid velocity.  
%                         These should be input as matrices. 
%          Fx,Fy          The x and y components of the forcing on the
%                         grid.  These are matrices the same size as U and
%                         V.
%          rho            The fluid density, which is a scalar.
%          mu             The diffusion coefficient, which is a scalar.
%          dx, dy, dt     The spatial step-size. The value for a FULL time-step.
%          Nx, Ny         The number of grid points along each side of the
%                         domain.
%          Lx, Ly         The length of the domain.
%          IndX, IndY     Matrix containing indices in x- and y-direction.
%                         Used as the wave number in FFT solution.
%
% OUTPUTS: U, V   The updated x and y components of the fluid
%                         velocity at a full time-step.  
%                         These will be in matrix format.
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%


% Construct FFT Operator
A_hat = 1 + 2*mu*dt/rho*( (sin(pi*IndX/Nx)/dx).^2 + (sin(pi*IndY/Ny)/dy).^2 );

% Compute first and second derivatives (centred).
Uxx = Centered_Dxx(U,dx);
Uyy = Centered_Dyy(U,dy);
Vxx = Centered_Dxx(V,dx);
Vyy = Centered_Dyy(V,dy);


% Construct right hand side in linear system
rhs_u = dt/rho*(  Fx + mu*(Uxx+Uyy) ...
        ) + U;

rhs_v = dt/rho*(  Fy + mu*(Vxx+Vyy) ...
        ) + V;

% Perform FFT
rhs_u_hat = fft2(rhs_u);
rhs_v_hat = fft2(rhs_v);  

% Calculate Fluid Pressure
p_hat = -(1.0i/dx*sin(2*pi*IndX/Nx).*rhs_u_hat + 1.0i/dy*sin(2*pi*IndY/Ny).*rhs_v_hat)./...
            ( dt/rho*((sin(2*pi*IndX/Nx)/dx).^2 +(sin(2*pi*IndY/Ny)/dy).^2 ));
        
% Zero out modes.
p_hat(1,1) = 0;
p_hat(1,Nx/2+1) = 0;
p_hat(Ny/2+1,Nx/2+1) = 0;  
p_hat(Ny/2+1,1) = 0;

% Calculate Fluid Velocity
u_hat = (rhs_u_hat - 1.0i*dt/(rho*dx)*sin(2*pi*IndX/Nx).*p_hat)./A_hat;
v_hat = (rhs_v_hat - 1.0i*dt/(rho*dy)*sin(2*pi*IndY/Ny).*p_hat)./A_hat;

% 2d inverse of Fourier Coefficients
U = real(ifft2(u_hat));
V = real(ifft2(v_hat));