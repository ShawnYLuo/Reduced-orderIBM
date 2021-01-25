function U = U0_shear(Y,ub)
y0 = Y(1,1)+Y(end,1);
U = 4*ub/y0/y0*Y.*(y0-Y);