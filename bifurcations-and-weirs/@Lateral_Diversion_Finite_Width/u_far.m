% Thu 16 Jan 20:06:07 +08 2020
function [uf] = u_far(obj,x,y)
	alpha = obj.alpha;
	gamma = obj.gamma;

	% radius
	if (~issym(x))
		r  =  hypot(x,y);
	else
		r = sqrt(x.^2 + y.^2);
	end
	uf =  obj.uin - alpha.*gamma.*x./(2*r).*coth(pi*gamma.*r./2);
end

