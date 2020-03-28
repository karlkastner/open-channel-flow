% Fri  9 Mar 17:56:18 CET 2018
%% along-channel surface elevation for (critical) lateral outflow over a side-weir
% assume zs1 = zs2 constant bed level
function [x, zs] = surface_elevation(obj,x,formflag);
	if (nargin() < 2 || isempty(x))
		x = linspace(0,obj.width.side)';
	end

	g     = obj.g;
	alpha = obj.alpha;

	Qu = obj.Q.upstream;
	hu = obj.h.upstream;
	wu = obj.width.upstream;
	Qs = obj.Q.side;
	ws = obj.width.side;
	%Au = hu*wu;

	switch (obj.method)
	case {'exact'}
		zs = hu/2*log(1 + (Qs^2*x.^2 - 2*Qs*Qu*ws*x) / ...
                                         (Qu^2*ws^2 - g*hu^3*ws^2*wu^2) ...
                                    );
	case {'approximal'}
		zs =  hu/2*(Qs^2*x.^2 - 2*Qs*Qu*ws*x) ...
                                / (Qu^2*ws^2 - g*hu^3*ws^2*wu^2);
	case {'numerical'}
		X = [0,ws];
		[x_,zs_] = ode45(@(x,z) obj.dzs_dx(x,z), X, 0);
		zs = interp1(x_,zs_,x);
	otherwise
		error('here');
	end
	zs = zs-zs(end);
end

