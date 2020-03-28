% Mon 10 Feb 17:46:48 +08 2020
% linearly varying velocities
% u1 : without inflow vel!
function [u,v] = velocity_linear(obj,x,y)
	gamma       = obj.gamma;
	alpha_dummy = 1;

%	[u0,v0] = obj.velocity(x,y,alpha_dummy,gamma);

	u1 = ( y.*v0 + x.*u0 -1);
	v1 = (-x.*v0 + y.*u0);

%	u1 = a.*u0 + b.*( y.*v0 + x.*u0-1);
%	v1 = a.*v0 + b.*(-x.*v0 + y.*u0);

	m = obj.m;
	dw = 1./m;

	iA = (spdiags(ones(m,1)*[0.5,0.5],0:1,5,m+1));
	iB = 1/dw*(spdiags(ones(5,1)*[-1,1],0:1,m,m+1));

	u = u0*IA*cp + u1*IB*cp;
	v = v0*IA*cp + v1*IB*cp;
end

