% Wed  8 Jan 18:20:09 +08 2020
function C = z02chezy(z0,h)
	g = Constant.gravity;
	kappa = Constant.Karman;
	C  = sqrt(g)/kappa*(log(h./z0) - 1);
	C(h<0) = NaN;
end

