% Thu  9 Aug 16:41:26 CEST 2018
%
%% load analytical solutions for potential flow field at a lateral diversion
%% with an infinitely wide main channel
%
% function [sol,obj] = lateral_outflow(obj,u0,Qs,ws,obj.shape)
%
function [sol,obj] = load_functions(obj)

	if (~isdeployed())
		what_   = what(class(obj));
		path_   = what_.path;
	else
		path_   = obj.deploypath;
	end
	load([path_,filesep(),obj.funfilename],'sol');

%	switch (length(varargin))
%	case {0}
%		% nothing to do, just return functions
%	case {1,2}
%		alpha = varargin{1};
%		if (length(varargin) > 1)
%			beta  = varargin{2};
%		else
%			beta = [];
%		end
		alpha = obj.alpha;
		beta  = obj.beta;

		% h = obj.h;
		obj.fun.phi = @(x,y)  sol.(obj.shape).nfun.phi(x,y,alpha,[]);
		obj.fun.u   = @(x,y)  sol.(obj.shape).nfun.u(x,y,alpha,[]);
		obj.fun.v   = @(x,y)  sol.(obj.shape).nfun.v(x,y,alpha,[]);
		obj.fun.ubed   = @(x,y)  sol.(obj.shape).nfun.ubed(x,y,alpha,beta);
		obj.fun.vbed   = @(x,y)  sol.(obj.shape).nfun.vbed(x,y,alpha,beta);
%		obj.fun.R   = @(x,y)  sol.(obj.shape).nfun.R(x,y,alpha,[]);
		obj.fun.J   = @(t,xy) sol.(obj.shape).nfun.J(xy(1),xy(2),alpha,beta);
		obj.fun.Jb  = @(t,xy) sol.(obj.shape).nfun.Jb(xy(1),xy(2),alpha,beta);
		obj.fun.y0  = @(x,y,alpha,beta)  sol.(obj.shape).nfun.y0([],[],alpha,[]);
%	otherwise
%		error('here');
%	end
end

