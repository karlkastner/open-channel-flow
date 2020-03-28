% Sun  9 Feb 17:13:17 +08 2020
% determine coefficients
function [cp, cs] = coefficients(obj)
	cp = obj.cp;
	if isempty(cp)
		v0 = -obj.alpha;
		u0 = obj.u0;
		n = obj.n;
		k = (-n:n);
		ns = obj.m;
		% segment mid:points
%  TODO the least squares oversampling fit does not go together with enforcing the mean of v

% TODO currently ill conditioned,
% match at more than one location per element and then due least squares?
% max value of du/dx ~ l^2 ~ pi^2 n^2, but potential velocity is zero at centres!


	m = ns;
	xl   = 1/(2*m)*((1-m):2:(m-1));
	xl = xl';
	dx = xl(2)-xl(1);
	x = xl;
if (0)
	x_ = [xl-dx/4;xl+dx/4];
else
	x_ = x;
end
	nn = length(x_);
%	x = -1/2+(1:ns)'/(ns+1)
		ns_ = ns/2;
	
		% non-zero frequencies
		%l = 2*pi*(1:ns);
		l = pi*(1:ns_);
		l = pi*(0:ns_-1);
		% matrix [pot,-a-fourier,-b-fourier]
		A = zeros(2*nn,2*ns);
		% rhs
		b = zeros(2*nn,1);
			
		% unweighted u value of vel at segment centres
		% unweighted u value 
		% Phi      =            (a*sin(l*x) + b*cos(l*x))*exp(-l*y)
		% u(x,0)   = dP/dx =    (   l*a*cos(l*x) - l*b*sin(l*x))
		% v(x,0)   = dP/dy =    (-l  *a*sin(l*x) - l  *b*cos(l*x))
		% du/dx    =            (-l^2*a*sin(l*x) - l^2*b*cos(l*x))
		% pot of u
		u_                         = evalk_('u',x_,0);
		A(     1:nn,      1:ns)    =  u_;
		%  - a of u
		A(     1:nn,   ns+1:ns+ns_)  =  -l.*cos(x_*l);
		%  - cos of u
		A(     1:nn, ns+ns_+1:2*ns)  =  +l.*sin(x_*l);
		% rhs u0
		b(     1:nn)               = -u0;
		% pot of v
		y = sqrt(eps);
		v_                         = evalk_('v',x_,y);
		A(  nn+1:2*nn,    1:ns)    =  v_;
		% - a of v
		A(  nn+1:2*nn,  ns+1:ns+ns_) = +l.*sin(x_*l);
		% - b of v
		A(  nn+1:2*nn,ns+ns_+1:2*ns) = +l.*cos(x_*l);
		% force mean velocity
		% TODO this is a quick hack (!), th choice of n_ is not arbitrary,
		% continuity is violated at one segment (!)
%sin(x_*l)
%mean(sin(x_*l))
%cos(x_*l)
%c = mean(cos(x_*l));
%c(1:2:end) = 0
%x_
%mean(cos(pi*x_))
%mean(cos(2*pi*x_))
%mean(cos(3*pi*x_))
%mean(cos(4*pi*x_))
%pause

if (0)
if (0)
		A(2*nn,1:ns) = 0;
		A(2*nn,ns+1:2*ns)   = mean(l.*sin(x_*l));
		A(2*nn,2*ns+1:3*ns) = mean(l.*cos(x_*l));
		b(2*nn) = -v0;
else
		n_ = round(1.5*nn);
		A(n_,ns+1:end) = 0;
		A(n_,1:ns) = mean(v_);
		b(n_) = v0;
end
end
		% rhs -v0
		b(  nn+1:2*nn)             = +v0;
if (0)
		% pot of du/dx
		du_dx_                     = evalk_('du_dx',x_,0);
		A(2*nn+1:3*nn,     1:ns)   = du_dx_;
		% - a of du/dx
		A(2*nn+1:3*nn,  ns+1:2*ns) = (l.*l).*sin(x_*l);
		% - b of du/dx
		A(2*nn+1:3*nn,ns+ns_+1:2*ns) = (l.*l).*cos(x_*l);
		% rhs derivative match
		b(2*nn+1:3*nn)             = 0;
end
% force odd cos-coefficients to zero
if (0)
	l_ = 1:2:m;
	for idx=1:length(l_)
		A(end-idx+1,:) = 0;
		b(end-idx+1,:) = 0;
		A(end-idx+1,2*ns+l_(idx)) = 1;
	end
end

	x
	AA = A(1:ns,1:ns);
	cond(AA)
	AA = A(ns+1:end,ns+1:end)
	cond(AA)
	AA = A(ns+1:end-1,ns+1:end-1)
	cond(AA)
	AA = A(ns+1:end-2,ns+1:end-2)
	cond(AA)
pause
if (0)
x_
	AA = A(ns+1:2*ns,ns+1:2*ns)
	[Q,R] = qr(AA)
%figure()
%plot(AA)
	AA=AA*diag(1./colnorm(AA))
	colnorm(AA)
pause

	AA = A(2*ns+1:3*ns,2*ns+1:3*ns)
	[Q,R] = qr(AA)
	AA = A(1*ns+1:3*ns,1*ns+1:3*ns)
	[Q,R] = qr(AA)
	cond(A(ns+1:end,ns+1:end))
	cond(A(ns+1:2*ns,ns+1:2*ns))
	cond(A(2*ns+1:end,2*ns+1:end))
pause
end
if (0)
	% as this is potential flow, fix the value at 0
	cdx = ns+ns/2;%-round(0.5*ns);
	A(cdx,:) = 0;
	b(cdx,:) = 0;
	A(cdx,  ns+1:ns+ns/2) = sin(0.*l);
	A(cdx,ns+ns/2+1:2*ns) = cos(0.*l);
end
%wx
%sin(x*l)
%cos(x*l)
%pause
%l.*l	
%figure(10)
%clf
%imagesc((A))
	
		% solve for coefficients
		cond(A)
		%c = A \ b;
		c = A \ b;
		%[A*c,b]
		%c = pinv(A)*b;
%c
%pause
%size(c)
%size(A)
%size(b)
%pause
		%c = pinv(A) * b;
%[A*c, A(:,1:ns)*c(1:ns),b-A(:,1+ns:end)*c(1+ns:end)]
%[b,A*c, A(:,1:ns)*c(1:ns),A(:,1+ns:2*ns)*c(1+ns:2*ns),A(:,2*ns+1:3*ns)*c(2*ns+1:3*ns)]
%pause
		% write back
		cp = c(1:ns);	
		cs = c(ns+1:end);
		obj.cp = cp;
		obj.cs = cs;
	end % if not isempty cp

function val = evalk_(fname,xi,y)
	alpha_ = obj.alpha;
	xx  = (xi-x');
	val1 = sum(obj.fun.k.(fname)(flat(xx)*ns, y*ns ...
			        , alpha_,[],obj.gamma ...
				, k ...
			       ),2);
	val = reshape(val1,size(xx));
%	val = sum(val);
end % evalk_

end % coefficients

