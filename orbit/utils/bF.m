function [Fx, Fy] = bF(sigma, x,y,x0,y0)
% The membrane position is X=(x,y). b0 and a0 are cosines and sines of the 
% angles of adjacent edges of the initial configuration.

[a, b] = tri_param(x,y);
[a0, b0] = tri_param(x0,y0);
par = a0.*b./a-b0;


[fx1,fy1] = patcos(x,y);

Fx = -par.*fx1;
Fy = -par.*fy1;


Fx = sigma*Fx;
Fy = sigma*Fy;

