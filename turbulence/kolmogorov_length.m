function eta = kolmogorov_length(eps)
	nu  = Physics.viscosity_kinematic_water();
	eta = (nu.^3./eps).^(1/4);
end

