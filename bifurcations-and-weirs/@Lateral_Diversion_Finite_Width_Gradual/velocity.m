% Thu 16 Jan 20:06:07 +08 2020
function [u,v] = velocity(obj,x,y,alpha,gamma)
	%[u,v] = lateral_outflow_finite_width_(x,y,obj.alpha,obj.gamma,obj.n);

	if (nargin()<4)
		alpha = obj.alpha;
		gamma = obj.gamma;
	end

	u0 = obj.u0;

	% correct for outflow
	v  = 0;
	uf = 0;
	vf = 0;

	if (issym(x))
		syms k 
		u = 0;
	else
		k = -n:n;
		u  = u0;
	end

	for k=k
		% shift origin
		y_ = y - 2*k./gamma;
		u  = u - alpha./(2*pi)*log(    (1 + 4*x.^2 + 4*y_.^2 + 4*x) ...
                                            ./ (1 + 4*x.^2 + 4*y_.^2 - 4*x) );
		v  = v - alpha./pi*atan2( 4*y_, 4*x.^2 + 4*y_.^2 - 1);
	end
end

