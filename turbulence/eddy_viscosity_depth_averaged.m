% Mon 21 Oct 15:01:18 PST 2019
function nu_t = eddy_viscosity_depth_averaged(h,us)
	% TODO make constant a class that treats symbols
	if (issym(h))
		syms kappa	
	else
		kappa = 0.41;
	end
	nu_t = 1/6*kappa*us.*h;
end

