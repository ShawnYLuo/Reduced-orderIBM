function myPlotMembrane(X, Y, U, V,chiX,chiY,figNum)
% Plots the membrane
%
% Assumption: Assuming chiX and chiY are column vectors
% Assumption: Assuming chiX(i+1)-chiX(i) < .5 and chiY(i+1)-chiY(i) < .5, for all points that don't cross the boundary
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

if nargin==7
    figure(figNum)
else
    figure;
end
Lx = X(1,end);
Ly = Y(end,1);
quiver(X,Y,U,V);
axis([0 Lx 0 Ly]);
daspect([1 1 1])
hold on
plot([chiX; chiX(1)],[chiY; chiY(1)],'r');
drawnow
hold off

set(gca,'Box','on');
