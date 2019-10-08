% Wed  4 Jan 12:21:22 CET 2017
% after yang 2006
%
%% vertical profile of the streamwise velocity in non-uniform flow
function [u, a, b, uh] = vertical_velocity_profile(us,z0,h,dh_dx,z,flag)
	if (nargin()<6)
		flag = false;
	end
	kappa = Constant.KAPPA;

	% xu is not exacly known, approximation by log law
	uh = (us/kappa).*log(h/z0);
	a = 1/kappa*dh_dx.*uh./us;
	b = -(uh./us).^2.*dh_dx;
	xi = z./h;

	xi0 = z0/h;
%	xi0=z0;
	switch (flag)
	case {0} % parabolic eddy viscosity, numerically integrated
		dxi = xi(2)-xi(1);
%		dxi = z(2)-z(1);
		I = cumsum(((1-xi).^a)./xi)*dxi;
		%xi0 = z0;
		u = us/kappa*( (1-xi).^-a .* ( I + b/a - log(xi0) ) - b/a );

	case {1} % series expansion
		I = log(z) - a*xi + a*(a-1)/4*xi.^2;
		I = I-I(1);
		u = us/kappa*( (1-xi).^-a .* (I + b/a - log(xi0)) - b/a );
		%u = us/kappa*( (1-xi).^-a .* (log(z) - a*xi + a*(a-1)/4*xi.^2 + b/a - log(z0)) - b/a );
	case {2} % 
		% eq 34 in yang
		u = us/kappa*(  log(xi) - b*log(1-xi) ...
				-4*b/pi*log(1./cos(pi*xi/2) + tan(pi*xi/2)) ...
				-log(z0/h) );
	end

%	ximax = max(1 + u.*v/us^2)/(1-b);

%	clf
%	I_ = (log(xi)-a*xi+a*(a-1)/4*xi.^2)
%	A=[I_-I_(1),(I)];
%	o=median(A(:,1)-A(:,2))
%	plot([A(:,1) A(:,2)+o])
%	median(A(:,1)./A(:,2))
%	pause

end

