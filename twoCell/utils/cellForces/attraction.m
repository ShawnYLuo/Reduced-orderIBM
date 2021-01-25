function [f1x,f1y,f2x,f2y] = attraction(x1,y1,x2,y2,s,ld)
%attraction computes the  force between two nodes. Attraction force
%generates when the distance between two nodes,d, in [ld,inf], and is 
%given by                 
%                           f=sd^2.
% Inputs: (x1,y1)is the coordinate of node 1
%         (x2,y2)is the coordinate of node 2
%         s is the force spring coefficient
%         ld is the separation distance
%Outputs: (f1x,f1y) is the force on node (x1,y1)
%         (f1x,f1y) is the force on node (x1,y1)

dx=x1-x2;
dy=y1-y2;
d=norm([dx;dy]);
%f = s*((d>ld)*(d-ld)-1.2*(d<ld));
%f = s*(d>ld)*(d/(d-ld));
flag = abs(d-ld)>ld/2;
f = s*(flag*abs(d-ld)+ld/2)*(d>ld);
f1x = -dx/d*f;
f2x = -f1x;
f1y = -dy/d*f;
f2y = -f1y;