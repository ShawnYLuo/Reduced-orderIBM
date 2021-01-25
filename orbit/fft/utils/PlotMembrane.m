function PlotMembrane(X, Y, U, V, chiX, chiY, chiX0, chiY0, figNum)
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

Lx = X(1,end)+X(1,1);
Ly = Y(end,1)+Y(1,1);



clf;
axis([0 Lx 0 Ly]);
xlabel('x'); ylabel('y');

hold all;

quiver(X,Y,U,V);

plot([chiX;chiX(1)],[chiY;chiY(1)],'r');
plot([chiX0;chiX0(1)],[chiY0;chiY0(1)],'k');


axis square;

drawnow;
hold off;

set(gca,'Box','on');
