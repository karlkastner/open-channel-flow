% 2022-04-06 18:25:33.884537230 +0200
function S = laminar_to_turbulent_transition_slope(D,n)
	g  = 9.81;
	nu = Constant.viscosity_kinematic_water;
	% 1/n*D.^(2/3)*sqrt(S) = 1/3*g/nu*S.*D.^2;
	S = (3/n*D.^(2/3-2)*nu/g).^2;
end
