function myPlotMembrane(X, Y, U, V,chiX1,chiY1,chiX2,chiY2)
% Plots the membrane
%
% Assumption: Assuming chiX and chiY are column vectors
% Assumption: Assuming chiX(i+1)-chiX(i) < .5 and chiY(i+1)-chiY(i) < .5, for all points that don't cross the boundary
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%


figure(1)
Lx = X(1,end);
Ly = Y(end,1);
U = (U(1:end-1,:)+U(2:end,:))/2;
V = (V(:,1:end-1)+V(:,2:end))/2;
quiver(X,Y,U,V);
axis([0 Lx 0 Ly]);
daspect([1 1 1])
hold on
plot([chiX1; chiX1(1)],[chiY1; chiY1(1)],'r');
plot([chiX2; chiX2(1)],[chiY2; chiY2(1)],'r');
drawnow
hold off

set(gca,'Box','on');

