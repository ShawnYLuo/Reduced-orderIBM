function [f1x,f1y,f2x,f2y] = repulsion(x1,y1,x2,y2,a,cd)
%repulsion computes the repulsion force between two nodes. Repulsion force
%generates when the distance between two nodes,d, is less than the 
%separation distance, cd, and is given by                 
%                           f=a*d*(cd^2-d^2).
% Inputs: (x1,y1)is the coordinate of node 1
%         (x2,y2)is the coordinate of node 2
%         a is the force coefficient
%         cd is the separation distance
%Outputs: (f1x,f1y) is the force on node (x1,y1)
%         (f1x,f1y) is the force on node (x1,y1)

dx=x1-x2;
dy=y1-y2;
d=norm([dx;dy]);
%f = a*d*(cd^2-d^2);
%f = a*((d<cd)*(cd-d)-1.2*(d>cd));
flag = abs(d-cd)>cd/2;
f = a*(flag*abs(d-cd)+cd/2)*(d<cd);
f1x = dx/d*f;
f2x = -f1x;
f1y = dy/d*f;
f2y = -f1y;



