% Wed 28 Feb 12:04:05 CET 2018
%
%% across channel derivative of flow depth for a meandering river
% f is not a friction factor, but the shape of the near bank velocity profile
% TODO use shiono knight to compute f
% waning: d_mm was d_m
function dh_dr = dh_dr(obj,r,h,d_mm,f,u_bar)
	g      = obj.g;
	Rc     = obj.Rc;
	C      = obj.C;
%	u_bar  = obj.u_bar;
	us_bar = sqrt(g)/C*u_bar;

	rhos = Constant.density.quartz;
	rhow = Constant.density.water;

	D = 1e-3*d_mm(r);
	f = f(r);

	tau_c  = critical_shear_stress(1e3*D,obj.T_C);

	sqrt_shields = us_bar/sqrt((rhos-rhow)/rhow*g*D);
	Fr = f.*r./Rc;
	% ikeda 1987, 23
	dh_dr  = Fr.^2.*h./Rc*sqrt_shields ...
		 .* (-1.941 + 1.226*u_bar/us_bar) ./ (5.382*sqrt(tau_c));
end % dh_dr

