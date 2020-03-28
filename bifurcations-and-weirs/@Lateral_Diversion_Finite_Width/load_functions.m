% Mon 20 Jan 11:35:02 +08 2020
function fun = load_functions(obj)
	sol = obj.derive();

% todo the load from file should be here, derive should not save, or?

	obj.fun   = sol.fun;
	obj.fun.J = @(t,xy) obj.J(xy(1:end/2),xy(end/2+1:end));
	obj.fun.Jb = @(t,xy) obj.Jb(xy(1:end/2),xy(end/2+1:end));
	% obj.fun.phi = @(x,y)  sol.(shape).nfun.phi(x,y,alpha,[]);
%	obj.fun = sol.fun;
%	obj.fun.u      = @(x,y)  sol.fun.u(x,y,obj.alpha,[],obj.gamma,obj.n,x,y);
%	obj.fun.v      = @(x,y)  sol.fun.v(x,y,obj.alpha,[],obj.gamma,obj.n,x,y);
%	obj.fun.ubed   = @(x,y)  sol.fun.ubed(x,y,obj.alpha,obj.beta,obj.gamma,obj.n,x,y);
%	obj.fun.vbed   = @(x,y)  sol.fun.vbed(x,y,obj.alpha,obj.beta,obj.gamma.obj.n,x,y);

%	obj.fun.R   = @(x,y)  sol.(shape).nfun.R(x,y,alpha,[]);
%	obj.fun.J   = @(t,xy) sol.(shape).nfun.J(xy(1),xy(2),alpha,beta);
%	obj.fun.Jb  = @(t,xy) sol.(shape).nfun.Jb(xy(1),xy(2),alpha,beta);
%	obj.fun.y0  = @(x,y,alpha,beta)  sol.(shape).nfun.y0([],[],alpha,[]);
end

