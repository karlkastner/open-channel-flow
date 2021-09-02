% 38 in merckelbach 2006
function eps = eddy_diffusivity(z,us,h)
	kappa = 0.41;
	eps = us.^3./(h*kappa).*(h-z)./z;
end


