% Thu  9 Aug 16:41:26 CEST 2018
%%
%% analytical potential flow solution to lateral outflow from an infinitely
%% wide channel
function [u,v,phi] = lateral_outflow_analytic(obj,u0,v0,ws,shape,x,y)
	w = what('Potential_Flow');
	load([w.path,filesep(),'lateral-outflow-functions.mat']);

	if (0 == ws)
		% dirac delta
		u   = s.delta.fun.u(x,y);
		v   = s.delta.fun.v(x,y);
		phi = s.delta.fun.phi(x,y);
	else
		% velocity distribution
		u = s.(shape).fun.u(ws,x,y);
		v = s.(shape).fun.v(ws,x,y);
		phi = s.(shape).fun.phi(ws,x,y);
	end
	% scale up and add offset
	u = u0 + v0*u;
	v = v0*v;
end

