function myPlotMembrane(X, Y, U, V,chiX,chiY,chiXt,chiYt,figNum)
% Plots the membrane
%
% Assumption: Assuming chiX and chiY are column vectors
% Assumption: Assuming chiX(i+1)-chiX(i) < .5 and chiY(i+1)-chiY(i) < .5, for all points that don't cross the boundary
%
% Authors: Jeffrey Wiens and Brittany Froese, Copyright 2011-2012
%

if nargin==9
    figure(figNum)
else
    figure;
end
Lx = X(1,end);
Ly = Y(end,1);
U = (U(1:end-1,:)+U(2:end,:))/2;
V = (V(:,1:end-1)+V(:,2:end))/2;
quiver(X,Y,U,V);
axis([0 Lx 0 Ly]);
daspect([1 1 1])
hold on
h(1) = plot([chiX; chiX(1)],[chiY; chiY(1)]);
h(2) = plot([chiXt; chiXt(1)],[chiYt; chiYt(1)],'k');
drawnow
hold off
%legend(h([1,2]),{'num','true'},'Orientation','horizontal')
set(gca,'Box','on');
