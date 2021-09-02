% 2021-08-21 15:27:14.979176511 +0200
function Q = flow_rate(dp,A,L)
	mu = Constant.viscosity_dynamic_water;
	Q = dp*A^2/(8*pi*mu*L);
end
