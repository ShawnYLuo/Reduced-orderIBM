function M = dataInterp2D(X,Y,V,chiXl,chiYl,chiXr,chiYr)
% Given some sample data of 2D function V(X,Y) in meshgrid
% format, this function computes the matrix M whose (i,j) entry M(i,j) =
% V(chiX(i)-chiX(j),chiY(i)-chiY(j)) using built-in function interp2.
%
% Note: If query points (chiX(i)-chiX(j),chiY(i)-chiY(j)) exceed the
%       domain given by data (X,Y), interp2 returns NaN.
%
% Inputs: (X,Y,V): matrices of sample data in meshgrid format.
%         chiX,chiY: vectors of same size representing query points.
% Outputs: M: square matrix whose size is the same as the size of chiX and
%             chiY.
%
% Yushuang Luo 6/30/2020
%


% covert chiX and chiY into column vectors
chiXl = chiXl(:);
chiYl = chiYl(:);
chiXr = chiXr(:);
chiYr = chiYr(:);

% construct matrices Xq and Yq. Xq(i,j) = chiX(i)=chiX(j), similar for Yq.
Xq = chiXl-chiXr';
Yq = chiYl-chiYr';

% call MATLAB built-in function "interp2" to compute M
M11 = interp2(X,Y,V{1,1},Xq,Yq); % linear interpolation by default
M12 = interp2(X,Y,V{1,2},Xq,Yq); % linear interpolation by default
M21 = interp2(X,Y,V{2,1},Xq,Yq); % linear interpolation by default
M22 = interp2(X,Y,V{2,2},Xq,Yq); % linear interpolation by default
M = [M11 M12;M21 M22];