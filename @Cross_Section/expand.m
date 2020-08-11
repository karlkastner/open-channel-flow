% 2020-04-16 22:25:50.622465667 +0800
function zbi = expand(obj,xi)
	switch (obj.mode)
	case {'piecewise-linear'}
		cy  = obj.c;
		n = length(cy);
		%cx  = 0.5*(0:n-1)/(n-1);
		cx  = (0:n-1)/(n-1);
		zbi = interp1_piecewise_linear(cx,cy,xi);
	case {'polynomial'}
		c = obj.c;
	switch (obj.order)
		% TODO use vander
		case {'quadratic'}
			zbi = c(2)*xi.^2 + c(1);
		case {'quartic'}
			zbi = c(3)*xi.^4 + c(2)*xi.^2 + c(1);
		case {'sextic'}
			zbi = c(4).*xi.^6 + c(3)*xi.^4 + c(2)*xi.^2 + c(1);
		otherwise
		if (order<2)
			zbi    = zb_min + xi_.*(zb_max-zb_min);
		else
			zbi = xi_.^2*c(3) + xi_*c(2) + c(1);
		end
	end % switch order
	end % switch method
end

