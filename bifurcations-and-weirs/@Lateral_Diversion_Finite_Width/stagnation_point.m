% Fri 24 Jan 09:57:42 +08 2020
function x = stagnation_point(obj,x0,opt)
	if (nargin()<2)
		x0 = 0.5+sqrt(eps);
	end
	LB = 0.5+sqrt(eps);
	if (nargin()<3)
		opt = optimset('MaxIter',32);
	end
	x  = fzero_newton(@(x) obj.u(x,0), ...
			  @(x) obj.evalk('du_dx',x,0), ...
			  x0, opt, LB);

%	x_     = fzero_bisect(@(x) obj.u(x,0),(0.5+eps)*ones(size(x)),0.6*ones(size(x)),opt);
%	[x,x_]
%pause
%%	fdx    = isnan(x);
%	x(fdx) = x_(fdx);
%	x0(jdx,idx) = fzero(@(x) lateral_outflow_finite_width2(x,y,a(jdx),w0(idx),n),x0_);
end

