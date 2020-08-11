% 2020-04-16 22:14:24.958432483 +0800
function obj = init(obj)
	zb = obj.zb.c.param;
	switch (obj.mode)
	case {'piecewise-linear'}
		c = zb;
	case {'polynomial'}
	switch (obj.order)
	case {'quadratic'} % without linear term
		order_ = [1,0,1];
		dof = 2;
	case {'quartic'} % without linear or cubic term
		order_ = [1,0,1,0,1];
		dof = 3;
	case {'sextic'} % without odd terms
		order_ = [1,0,1,0,1,0,1];
		dof = 4;
	otherwise
		order_ = order;
		dof = length(order)+1;
	end
		% points where the bed level is given
		x0  = linspace(0,0.5,dof)';
		A   = vander_1d(x0,order_);
		% coefficients of the polynomial
		c   = A \ cvec(zb);
	end
	obj.c = c;
end

