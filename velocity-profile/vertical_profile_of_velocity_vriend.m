% 2018-08-30 11:00:47.531470407 +0200
%% vertical profile of the streamwise velocity, method of de vriend
function [ut,v,fs,fm,zt] = vertical_profile_of_velocity_vriend(ubar,d,R,R0,C,z,dubar_dphi)
	if (~issym(z))
		g     = Constant.gravity;
		kappa = Constant.Karman;
	else
		syms kappa g
		assume(kappa > 0)
		assume(g > 0)
	end
	s =  0;
	b = -1;

	z = z/d;
	
	% TODO this remains unclear how dudphi is calculated
	%dubar_dphi = 0;

	delta = d/R0

	% 41
	fm = 1 + sqrt(g)/(kappa*C) + sqrt(g)/(kappa*C)*log(1+z);

	% 43b f.
	zt = exp(-1 - kappa*C/sqrt(g));

	% 43a 
%	F_1 = int(log(1+z_)/z_,z_,-1+zt,z);
	F_1 = @(z) dilog(-z) + log(-z).*log(z + 1);
	F_1 = F_1(z) - F_1(-1+zt);

%	F_2 = int(log(1+z_^2)/z_,z_,-1+zt,z);
	F_2 = @(z) -polylog(2, -z.^2)/2;
	F_2 = F_2(z) - F_2(-1+zt);

	% 42
	fs = 2*F_1 + sqrt(g)/(kappa*C)*F_2 - 2*(1-sqrt(g)/(kappa*C)).*fm;

	% eq 44
	ut = ubar*fm - delta*(s-b)*1/kappa^2*dubar_dphi*fs;

	% eq 45
	% eq 48
	%kappa = 0.41;
	% v = -delta (s-b) 1/kappa^2 \bar u'/r fs
	r  =  R/R0;
	v  = -delta*(s-b)*1/kappa^2*ubar/r*fs;

if (0)
	% the noramlization in the paper is incomplete and unclear, but more
	% elaborated on in the report, c.f. eq. 109 therein:
	eps  = delta % = d/R
	chib = b; % = zb/d
	chi  = s; % z/d  % = s?
	%s    = z/d
	r    = R/R0
	v    = -eps*(chi-chib)*1/kappa^2*umag/r*fs
end
end

