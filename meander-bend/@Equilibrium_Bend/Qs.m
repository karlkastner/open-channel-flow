% 2020-12-06 17:06:20.140098298 +0100
function Qs = Qs(obj,h)
	if (nargin()<2)
		h = obj.h;
	end
	u  = obj.u(h);
	% TODO no magic numbers
	% TODO apply f
	% TODO apply 1/theta slope
	p  = 5;
	D  = 1e-3*obj.D;	
	Qs = mean((u.^p)./obj.D)*obj.width;
end

