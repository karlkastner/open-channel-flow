function tau = kolmogorov_time(eps)
	nu  = Physics.viscosity_kinematic_water();
	tau = (nu./eps).^(1/2);
end

