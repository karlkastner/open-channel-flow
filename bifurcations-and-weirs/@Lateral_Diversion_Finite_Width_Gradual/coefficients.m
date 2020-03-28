% Sun  9 Feb 17:13:17 +08 2020
% determine coefficients
function [cp, cs] = coefficients(obj)
	mode = obj.mode;
	lmode = obj.lmode;
	ox = 0.5;

	cp = obj.cp;
	if isempty(cp)
		v00 = -obj.alpha;
		u00 = obj.u0;
		ns = obj.m;

		m = ns;

		% matching locations
		xmatch  = -1/2+(1:m)'/(m+1);

		[xp, dw, np] = obj.xp;

		nn = length(xmatch);
		ns_ = ns;
	
		% non-zero frequencies
		if (lmode)
			% force odd coefficients to zero
			l = pi*(2:2:(ns_));
	%		w_ = (1-l/(pi*(2*ns_)))
	if (obj.cflag)
	w_ = (1-l/(pi*(2*ns_)))
	else
	w_ = ones(size(l));
	end
		else
			% for more than 12 it becomes ill conditioned
			ns_ = min(10,ns_);
			l   = pi*(1:ns_);
			w_  = ones(size(l));
		end
		ns_ = length(l);
		% matrix [pot,-a-fourier,-b-fourier]
		A = zeros(2*nn,np+ns_);
		% rhs
		b = zeros(2*nn,1);
			
		% unweighted u value of vel at segment centres
		% unweighted u value 
		% it is exp(+y*l), not -yl here!
		% Phi      =            (a*sin(l*x))*exp(+l*y)
		% u(x,0)   = dP/dx =    (+l*a*cos(l*x))
		% v(x,0)   = dP/dy =    (+l*a*sin(l*x))
		% pot of u
		%ymatch = sqrt(eps)*ones(size(xmatch));
		ymatch = 10*(eps)*ones(size(xmatch));
		switch (mode)
		case {0}
			u_                         = obj.evalk_('u',xmatch,0.*ymatch);
			A(     1:nn,      1:np)    =  u_;
			v_                         = obj.evalk_('v',xmatch,ymatch);
			A(  nn+1:2*nn,    1:np)    =  v_;
		case {1}
			u0 = obj.evalk_('u',xmatch,ymatch);
			v0 = obj.evalk_('v',xmatch,ymatch);
		
			yp = 0;	
			dx = (xmatch-xp')./dw';
			dy = (ymatch-yp')./dw';
			u1 = ( dy.*v0 + dx.*u0 + 1/(pi)); % why +
			v1 = ( dx.*v0 - dy.*u0);
	%		u1 = ( ymatch.*v0 + xmatch.*u0 - 1);
	%		%v1 = (-xmatch.*v0 + ymatch.*u0);
	%		v1 = ( xmatch.*v0 - ymatch.*u0);

			IA = (spdiags(ones(np,1)*[0.5,0.5],0:1,np-1,np));
			%IB = 1/dw(1)*(spdiags(ones(np,1)*[-1,1],0:1,np-1,np));
			IB = (spdiags((1./dw.^0)*[-1,1],0:1,np-1,np));
			%[u0,v0,u1,v1]   = evalk1(xmatch,y0);
			A(1:nn,1:np)      = u0*IA + u1*IB;
			% pot of v
			A(nn+1:2*nn,1:np) = v0*IA + v1*IB;
		end
		%  - a of u
		A(     1:nn,   np+1:np+ns_)  = -w_.*l.*cos((xmatch)*l);
		% rhs u0
		b(     1:nn)                 = -u00;
		% - a of v
		A(  nn+1:2*nn,  np+1:np+ns_) = -w_.*l.*sin((xmatch)*l);
		% - b of v
		% rhs -v0
		b(  nn+1:2*nn)               = +v00;

		% solve for coefficients
%figure(10)
%clf
%imagesc(A)
%cond(A)
%pause
		c   = A \ b;
		%c   = pinv(A'*A) \ A'*b;
		res = A*c-b;
		fprintf('cond(A) %f\n',cond(A))
		fprintf('rms %f\n',rms(res))
		% write back
		cp     = c(1:np);
		cs     = c(np+1:end);
		obj.cp = cp;
		obj.cs = cs;
	end % if not isempty cp

end % coefficients

