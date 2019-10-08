% 2018-05-19 13:42:08.011870441 +0200
%
%% transverse gradient of the bed level of an equilibrium meander bend
%% for the case of uniform bed material
%	C      = 60;
%	d50_mm = 2*0.25;
% f : velocity profile (can account for near bank layer)
% for instantaneous profiles where u is a function : pass r=R
function dh_dr = dh_dr_uniform(obj,r,h,u_bar,d50_mm,f)
	g      = obj.g;
	Rc     = obj.Rc;
	C      = obj.C;
	rho_w  = Constant.density.water;
	f      = f(r);
	tau    = rho_w*g./C.^2*u_bar.^2;
	tau_c  = critical_shear_stress(d50_mm);

	% eq 8
	% velocity
	Fr = f*r./Rc;

	% ikeda, 1987, eq 21
	dh_dr =  Fr*h./r.*(-1.941 + 1.226*C/sqrt(g)) ...
                        ./(5.382*sqrt(tau_c./tau));
end % dh_dr

