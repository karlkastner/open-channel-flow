% 2018-02-09 10:04:33.847494082 +0100
% Wed 14 Mar 16:06:34 CET 2018
%% compute a streamline
% near bed streamlines
% function [t, xy] = streamline(obj, T, xy0, ufield, vfield, opt)
function [t, xy] = streamline(obj, T, xy0, ufield, vfield)
	if (nargin() < 6)
		% TODO choose on grid-size
		opt = struct();
	end
%	switch (ufield)
%	case {'u'}
%		% TODO this should be done in the inherited function 
%		if (isa(obj,'Potential_Flow_Analytic'))
%			J = obj.fun.J;
%			opt = odeset(opt,'jacobian',J);
%		end
%	case {'ubed'}
%		if (isa(obj,'Potential_Flow_Analytic'))
%			Jb = obj.fun.Jb;
%			opt = odeset(opt,'jacobian',Jb);
%		end
%	end
%	% TODO, choose on grid size
%	opt.InitialStep = 1e-3*T;
%	if (~isfield(opt,'MaxStep'))
%		opt = odeset(opt,'MaxStep',obj.streamlineopt.MaxStep);
%	end
%	if (nargin()>5)
%		opt=odeset(opt,'Events', stopevent);
%	end
% TODO, no magic numbers

	solver = obj.streamlineopt.solver;

	% TODO, nasty fix
	global xystop
	xystop = xy0;
	xy = [];
	t = [];
	while (~isempty(xystop))
		% allow bed flow leaving inner bend to reemerge at outer opposit bank
		xy0 = xystop;
		xystop = [];
		%solver = @ode15s;
		%solver = @ode45;
		[t_,xy_] = solver(@(t,xy) [obj.(ufield)(xy(1),xy(2));
					   obj.(vfield)(xy(1),xy(2))], ...
                                           T, xy0, obj.streamlineopt);
		%[t_,xy_] = ode23s(@(t,xy) phidot(xy), T, xy0,opt);
		t  = [t; NaN; t_];
		xy = [xy; NaN(1,2); xy_];
	end

%	function uv = phidot(xy)
%		uv = [obj.interp.(ufield)(xy(1),xy(2));
%		      obj.interp.(vfield)(xy(1),xy(2))];
%	end
end % streamline

