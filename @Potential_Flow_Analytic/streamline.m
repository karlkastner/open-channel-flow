% Sun 13 Jan 11:40:48 CET 2019
% Karl Kastner, Berlin
%
%% numerically follow path along streamline by integrating the velocity
function [t,xy] = streamline(obj, T, xy0, ufield, vfield, opt)
	if (nargin() < 6)
		opt = struct();
	end
	switch (ufield)
	case {'u'}
		opt = odeset(opt,'jacobian',obj.fun.J);
	case {'ubed'}
		if (~isempty(obj.fun.Jb))
			opt = odeset(opt,'jacobian',obj.fun.Jb);
		end
	otherwise
		error('here');
	end % ufield
	[t, xy] = streamline@Potential_Flow(obj, T, xy0, ufield, vfield, opt);
end % streamline
