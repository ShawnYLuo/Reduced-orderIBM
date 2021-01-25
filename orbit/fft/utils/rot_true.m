function [X,Y] = rot_true(X0,Y0,cx,cy,a)
trans = [cos(a) -sin(a);sin(a) cos(a)]*[X0(:)'-cx;Y0(:)'-cy]+[cx;cy];
X=trans(1,:)';
Y=trans(2,:)';