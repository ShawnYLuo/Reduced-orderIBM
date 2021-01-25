function [fx,fy] = elasticity(chiX,chiY,sigma,L,ds,Nb)
dChiX = (chiX([2:Nb,1])-chiX)/ds;
dChiY = (chiY([2:Nb,1])-chiY)/ds;

% We also want the magnitude of the derivative.
dChi = sqrt(dChiX.^2 + dChiY.^2);

% Now we compute the tension, which is also given at the points (s+1/2),
% and points in the direction of the tangent vector.
Tx = sigma * dChiX .* (1 - L ./ dChi);
Ty = sigma * dChiY .* (1 - L ./ dChi);

% The force density is computed at the points s by taking a narrow centred
% difference of the tension. This is equivalent to computing a backward
% difference of T_{s+1/2}.
% This assumes a simple linear spring model with resting length L.
% For the backward difference on the periodic domain, the value of s_0
% is equivalent to the value of s_{Nb}.
fx = (Tx - Tx([Nb,1:Nb-1])) / ds;
fy = (Ty - Ty([Nb,1:Nb-1])) / ds;