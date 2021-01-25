function U = U0_shear(Y,ub,Ly)
U = ub*(ones(size(Y))-2*Y/Ly);