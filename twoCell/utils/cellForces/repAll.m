function [f1x,f1y,f2x,f2y] = repAll(x1,y1,x2,y2,a,cd)
%repAll computes the repulsion force on two sets of nodes from each other. 
%Repulsion force generates when the distance between two nodes,d, is less 
%than the separation distance, cd, and is given by                 
%                           f=a/(d-cd).
% Inputs: x1 and y2 are vectors representing coordinates of some nodes
%         x2 and y2 are vectors representing coordinates of the other nodes
%         a is the force coefficient
%         cd is the separation distance
%Outputs: f1x and f1y are vectors of force on nodes (x1,y1)
%         f2x and f2y are vectors of force on nodes (x2,y2)


% point-to-point force
f1x = zeros(size(x1));
f1y = zeros(size(y1));
f2x = zeros(size(x2));
f2y = zeros(size(y2));
for i=1:length(x1)
    for j=1:length(x2)
        [f11,f12,f21,f22] = repulsion(x1(i),y1(i),x2(j),y2(j),a,cd);
        f1x(i)=f1x(i)+f11;
        f1y(i)=f1y(i)+f12;
        f2x(j)=f2x(j)+f21;
        f2y(j)=f2y(j)+f22;
    end
end

% same force
cx1 = mean(x1);
cy1 = mean(y1);
cx2 = mean(x2);
cy2 = mean(y2);
[f11,f12,f21,f22] = repulsion(cx1,cy1,cx2,cy2,a,cd);
f1x = f11*ones(size(x1));
f1y = f12*ones(size(y1));
f2x = f21*ones(size(x2));
f2y = f22*ones(size(y2));