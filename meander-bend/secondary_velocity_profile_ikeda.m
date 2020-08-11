% Sun  3 May 12:26:28 +08 2020
% c.f. ikeda, kikkawa 1985
% based on linear eddy viscosity profile
% parabolic eddy viscosity profile leads to infinite value near bed
function v = secondary_velocity_profile_ikeda(z,h,R,U,us,f,method)
	if (~issym(z))
		k = Constant.Karman;
	else
		syms k;
	end
	if (nargin()<7)
		method = 'exact';
	end
	switch (method)
	case {'exact'}
		xi = z/h;
		Fa = -15*(xi^2*log(xi) - 1/2*xi^2 + 15/54);
		Fb =  15/2*(xi^2*log(xi)^2 - xi^2*log(xi)+1/2*xi^2-19/54);
		v = f^2*U*h/R*1/k*(Fa - 1/k*us/U*Fb);
		% v = -(5*f^2*h*(30*U*k - 19*us + 27*us*xi^2 - 54*us*xi^2*log(xi) + 54*us*xi^2*log(xi)^2 - 54*U*k*xi^2 + 108*U*k*xi^2*log(xi)))/(36*R*k^2)
	case {'linear','linearized'}
		% this leads to the well known relations:
		% v|_0 = -12.0*U*f^2*h/R
		% and
		% z|v=0 = 0.47 h ~ 1/2 h
		%v = (25.0*U*f^2*z)/R - (12.0*U*f^2*h)/R;
		xi = z/h;
		% v = f^2*U*h/R*(25.0*xi - 12.0);
		%v = f^2*U/k^2*h/R*(-4.89*k + 1.3*us/U + (10.4*k - 3.6*us/U)*xi);
		v = f^2*U*h/R*(-11.9 + 7.73*us/U + (25.4 - 21.4*us/U)*xi);
%		v = (7.73*f^2*h*us)/R - (11.9*U*f^2*h)/R + (25.4*U*f^2*h*xi)/R - (21.4*f^2*h*us*xi)/R
	case {'cubic'}
		% maybe it is smarter to expand v at 0
		xi = z/h;
		% note : u0 cannot be ignored to get velocity near bed right
		%v= - (24.4*U*h*f^2*xi^3)/R + (25.4*U*h*f^2*xi^2)/R + (18.3*U*h*f^2*xi)/R - (11.7*U*h*f^2)/R;
		%v = f^2*U*h/R*(-24.4*xi^3 + 25.4*xi^2 + 18.3*xi - 11.7);
		v  = (20.1*f^2*h*us)/R - (11.7*U*f^2*h)/R + (18.3*U*f^2*h*xi)/R - (75.5*f^2*h*us*xi)/R + (25.4*U*f^2*h*xi^2)/R - (24.4*U*f^2*h*xi^3)/R + (67.8*f^2*h*us*xi^2)/R - (18.3*f^2*h*us*xi^3)/R;
	end
end

