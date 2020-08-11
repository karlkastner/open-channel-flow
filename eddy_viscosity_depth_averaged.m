% Mon 21 Oct 15:01:18 PST 2019
function nu_t = eddy_viscosity_depth_averaged(H,us)
	kappa = 0.41;
	nu_t = 1/6*kappa*us.*H;
end

