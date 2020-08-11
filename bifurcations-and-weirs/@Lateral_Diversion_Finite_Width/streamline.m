% Fri 24 Jan 18:32:21 +08 2020
function [t,xy] = streamline(obj,T,xy0,field,opt)
	if (nargin()<5)
		opt=odeset();
	end
	if (isempty(field))
		field = '';
	end
	switch (field)
	case {'bed'}
		opt = odeset(opt,'jacobian',obj.fun.Jb);
		fun = @uvbode;
	otherwise
		opt = odeset(opt,'jacobian',obj.fun.J);
		fun = @uvode;
	end
	% ode23s  : 53.694086
	% ode23tb : 59.595255
	% ode1ts  : 81.470467
	% ode23t  : 102.005731
	% ode45   : no completion in reasonable time

	%[t,xy] = ode23s(fun, T, xy0,opt);
	[t,xy] = ode23tb(fun, T, xy0,opt);

function uv   = uvode(t,xy)
	n     = length(xy)/2;
	x     = xy(1:n);
	y     = xy(n+1:2*n);
	[u,v] = obj.uv(x,y);
	uv = [u;v];
end

function uvb   = uvbode(t,xy)
	n     = length(xy)/2;
	x_     = xy(1:n);
	y_     = xy(n+1:2*n);
	[u,v] = obj.uvb(x_,y_);
	uvb = [u;v];
end

end
