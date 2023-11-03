
opt = optimset ( "MaxIter" , 50,"MaxFunEvals",50,'PlotFcns',@optimplotfval);
z0 = zeros(17,1);
z0(8,1)=1.0;
%chi=merit_sph(z0);
z = fminsearch(@merit_sph,z0,opt);