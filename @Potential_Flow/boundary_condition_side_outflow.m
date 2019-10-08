% Thu 17 May 14:31:02 CEST 2018
%% apply boundary conditions for side outflow
%% p*phi + (1-p)*d/db phi = rhs
%% y : along channel coordinate
function [p,rhs] = boundary_condition_side_outflow(obj,x,y,h,id,xlim,ylim,yc,ws,u0,qside,shapeflag)
	p   = zeros(size(x));
	rhs = zeros(size(x));
	h   = h(id);

	phi0 = u0*(ylim(2)-ylim(1));

	tol = sqrt(eps);

	% left, no flow across boundary
	fdx      = abs(x-xlim(1)) < tol;
	p(fdx)   = 0; % set phi' = v = 0
	rhs(fdx) = 0;

	% right, no flow across bank
	fdx      = abs(x-xlim(2)) < tol;
	p(fdx)   = 0; % set phi' = v = 0
	rhs(fdx) = 0;

	% right, flow into outlet
	fdx      = abs(x-xlim(2)) < tol & abs(y-yc) < ws/2;
	p(fdx)   = 0;
	switch(shapeflag)
	case {0,'const','rect'} % constatn
		rhs(fdx) = qside./h(fdx);
	case {1,'v'} % v
		dyrel       = max(0,2*(1-2*abs(y-yc)/ws));
		rhs(fdx) = qside./h(fdx).*dyrel(fdx);
	case {2} % ramp up
		dyrel    = 2*min(1,max(0,(y-yc+ws/2)/ws));
		rhs(fdx) = qside./h(fdx).*dyrel(fdx);
	case {3} % ramp down
		dyrel    = 2*(min(1,max(0,1+((yc-ws/2)-y)/ws)));
		rhs(fdx) = qside./h(fdx).*dyrel(fdx);
	case {'quad'}
		dyrel    = 30*((y-yc)/ws-1/2).^2.*((y-yc)/ws+1/2).^2;
		rhs(fdx) = qside./h(fdx).*dyrel(fdx);
	case {4,'cos'} % cos
		dyrel = 2*(y-yc)/ws;
		rhs(fdx) = qside./h(fdx).*(1+cos(pi*dyrel(fdx)));
	end
	
	% bottom, flow into main channel
	fdx      = abs(y-ylim(1)) < tol;
	p(fdx)   = 1; % set phi = int u
	rhs(fdx) = 1/2*phi0;

	% top, flow out of main channel
	fdx      = abs(y-ylim(2)) < tol;
	p(fdx)   = 1; % set phi = int u
	rhs(fdx) = -1/2*phi0;
end % boundary_condition_side_outflow

