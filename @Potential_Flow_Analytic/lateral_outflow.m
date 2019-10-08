% Thu  9 Aug 16:41:26 CEST 2018
%
%% potential flow solution to the case of lateral outflow from an infinitely
%% wide channel
%
% function [sol,obj] = lateral_outflow(obj,u0,Qs,ws,shape)
%
function [sol,obj] = lateral_outflow(obj,shape,varargin)
	w     = what(class(obj));
	path_ = [w.path,filesep(),obj.funfilename];
	load(path_,'sol');

	switch (length(varargin))
	case {0}
		% nothing to do, just return functions
	case {1,2}
		alpha = varargin{1};
		if (length(varargin) > 1)
			beta  = varargin{2};
		else
			beta = [];
		end
		% h = obj.h;
		obj.fun.phi = @(x,y)  sol.(shape).nfun.phi(x,y,alpha,[]);
		obj.fun.u   = @(x,y)  sol.(shape).nfun.u(x,y,alpha,[]);
		obj.fun.v   = @(x,y)  sol.(shape).nfun.v(x,y,alpha,[]);
		obj.fun.ubed   = @(x,y)  sol.(shape).nfun.ubed(x,y,alpha,beta);
		obj.fun.vbed   = @(x,y)  sol.(shape).nfun.vbed(x,y,alpha,beta);
%		obj.fun.R   = @(x,y)  sol.(shape).nfun.R(x,y,alpha,[]);
		obj.fun.J   = @(t,xy) sol.(shape).nfun.J(xy(1),xy(2),alpha,beta);
		obj.fun.Jb  = @(t,xy) sol.(shape).nfun.Jb(xy(1),xy(2),alpha,beta);
		obj.fun.y0  = @(x,y,alpha,beta)  sol.(shape).nfun.y0([],[],alpha,[]);
	case {3,4}
		u0 = varargin{1};
		Qs = varargin{2};
		ws = varargin{2};
		if (length(varargin) > 3)
			h  = varargin{3};
		else
			h = [];
		end
		obj.h([],[],h);
		obj.fun.phi = @(x,y)  sol.(shape).fun.phi(x,y,u0,Qs,h,ws);
		obj.fun.u   = @(x,y)  sol.(shape).fun.u(x,y,u0,Qs,h,ws);
		obj.fun.v   = @(x,y)  sol.(shape).fun.v(x,y,u0,Qs,h,ws);
		obj.fun.R   = @(x,y)  sol.(shape).fun.R(x,y,u0,Qs,h,ws);
		obj.fun.J   = @(t,xy) sol.(shape).fun.J(xy(1),xy(2),u0,Qs,h,ws);
		obj.fun.Jb  = @(t,xy) sol.(shape).fun.Jb(xy(1),xy(2),u0,Qs,h,ws);
	otherwise
		error('here');
	end
end

