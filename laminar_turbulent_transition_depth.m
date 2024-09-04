% 2022-04-06 18:31:00.793524048 +0200
function D = laminar_trubulent_transition_depth(S,n)
	g  = 9.81;
	nu = Constant.viscosity_kinematic_water;
	% 1/3*g/nu*S*D^2 = 1/n*D.^(2/3)*sqrt(S);
	 D = (3/g*nu*1./n*1./sqrt(S)).^(3/4);
end

