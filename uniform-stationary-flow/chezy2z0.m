% Wed  8 Jan 18:20:16 +08 2020
function z0 = chezy2z0(C,h)
	g     = Constant.gravity;
	kappa = Constant.Karman;
	z0    = h./exp((C*kappa/sqrt(g))+1);
end

