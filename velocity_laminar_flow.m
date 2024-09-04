% 2022-04-06 18:11:09.054664099 +0200
function u = velocity_laminar_flow(D,S)
	g  = 9.81;
	nu = Constant.viscosity_kinematic_water;
	u  = 1/3*g/nu*S.*D.^2;
end
 


