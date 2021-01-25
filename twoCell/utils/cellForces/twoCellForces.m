function [fx1,fy1,fx2,fy2] = twoCellForces...
    (sigma1,sigma2,s,ld,a,cd,x1,y1,x10,y10,x2,y2,x20,y20)
% oneCellForce computes the force density on two cells. The forces include
% bending force, repulsion and attraction.

[fx1b, fy1b] = bF(sigma1, x1,y1,x10,y10);
[fx2b, fy2b] = bF(sigma2, x2,y2,x20,y20);
[fx1a,fy1a,fx2a,fy2a] = attAll(x1,y1,x2,y2,s,ld);
[fx1r,fy1r,fx2r,fy2r] = repAll(x1,y1,x2,y2,a,cd);
fx1 = fx1b+fx1a+fx1r;
fy1 = fy1b+fy1a+fy1r;
fx2 = fx2b+fx2a+fx2r;
fy2 = fy2b+fy2a+fy2r;

% figure(2)
% quiver(x1,y1,fx1b,fy1b,'b');
% hold on
% quiver(x1,y1,fx1a,fy1a,'k');
% quiver(x1,y1,fx1r,fy1r,'r');
% quiver(x2,y2,fx2b,fy2b,'b');
% quiver(x2,y2,fx2a,fy2a,'k');
% quiver(x2,y2,fx2r,fy2r,'r');
% plot([x1; x1(1)],[y1; y1(1)],'m');
% plot([x2; x2(1)],[y2; y2(1)],'m');
% axis([min(x1)-0.1 max(x1)+0.1 2 4])
% daspect([1 1 1])
% drawnow
% hold off
