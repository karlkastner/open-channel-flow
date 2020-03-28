% 2013/04/14 23:45
% Karl Kästner, Berlin
%% predefined functions of channel geometries
% TODO move into swe class
function y = swe_ic(x, Xi, ictype, zbfun, wfun, cdfun, Q0, a, uflag)
	if (nargin()<9)
		uflag = false;
	end

%	g   = 9.81;
%	fps = 1/24;
%	nx = 2000;
%	t_end = 75;
%       T = 5;
%	L = 160;
%	L0 = 0.5*L;
%       dx = L/(nx-1);
%	x = dx*(-1:nx)' - L0; % + 2 boundary points
%	H0 = 1;
%	h0 = 0.5;
	
	L  = Xi(2)-Xi(1);
	x0 = 0.5*(Xi(2)-Xi(1));
	n  = length(x);

	switch (ictype)
	case {'flat'}
		h  = zeros(n,1);
	case {'hump'} % single gaussian wave
		% for hole just choose negative a
		h  = @(x) a*normpdf(3*(x-x0)/L);
		%h = H0*ones(nx+2,1) + h0*normpdf(x,0,0.025*L)/normpdf(0,0,0.025*L);
		%h = H0*ones(nx,1) - h0*normpdf(x,0,0.1*L);
	case {'wave-with-trough'}
		h = H0*ones(nx,1) + h0*(normpdf(x,0.45*L,0.1*L) - normpdf(x,0.55*L,0.1*L));
	case {'rectangle'}
		TODO
	case {'saw-tooth'}
		t = max(0, 1 - 0.5*0.06125*(x-x0)).*(x > x0);
		h = a*ones(nx+2,1) + a*t;
		%h = ones(nx+2,1) + 0.5*h0*(abs(x-L0) < L/10);
	case {'dambreak'}
		h  = @(x) a+a*(x<x0);
		% h = H0*ones(nx+2,1) + h0*(x < 0);
	case {'dambreak2'}
		h = a+a*(x<x0-1/6*L)+a*(x<x0+1/6*L);
	case {'wavepacked'}
		h = a*normpdf(0.5*w.*x).*cos(2*pi*x.*wg*0.5);
	case {'train'}
		h = a*0.125*sinc(x);
	case {'backwater'}
		g = 9.81;
		Cfun = @(x) sqrt(g./cdfun(x));
		zs0 = 0;
		[x_ h_ zs_] = backwater_curve(Q0,Cfun,wfun,zbfun,zs0,Xi);
		h = interp1(x_,zs_,x);
		% q = [h0__;-q0*ones(length(x),1)];
	otherwise 
		error('here')
	end
	
	% subtract bottom
	h = h - zbfun(x);

	% initial velocity
	switch (uflag)
		case {'right'} % initial velocity for single wave packet
			% make one characteristic become zero
			u = -sqrt(g)*(h-H0)/sqrt(H0);
			q = h.*u;
	%		u = - sqrt(g*h/h0);
		otherwise 
			% zero velocity, initial wave will split in travel in both directions
			% mean flow is offset
			q = -Q0./wfun(x);
%*ones(n,1);
	end
%	A = wfun(x).*h;
	y = [h; q];
end % swe_ic

