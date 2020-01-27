% Wed  8 Jan 18:20:09 +08 2020
function C = z02chezy(z0,h)
	g = 9.81;
	kappa = 0.41;
	C  = sqrt(g)/kappa*(log(h./z0) - 1);
end

