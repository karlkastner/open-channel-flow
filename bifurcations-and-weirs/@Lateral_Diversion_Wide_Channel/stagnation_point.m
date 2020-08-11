% Mon  6 Apr 16:03:35 +08 2020
function x0 = stagnation_point(obj,alpha)
	if (nargin()<2)
		alpha = obj.alpha;
	end
	x0 = obj.fun.y0([],[],alpha,[]);
end
