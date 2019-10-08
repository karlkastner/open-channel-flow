% Mo 6. Jul 08:48:59 CEST 2015
% Karl Kastner, Berlin
%
%% transverse (across channel) profile of the streamwise velocity
%% in a straight channel
%% numerical solution
%% the eps seems incorrect, use better stationary_1d_swe
%%
%% rho g h S - beta q^2 f / (8 h^2)   + d/dy(eps_t dq/dy) = 0
%% rho g h S - beta q^2 g / (C^2 h^2) + d/dy(eps_t dq/dy) = 0
function [u, q, delta] = lateral_division_method(N,h,C,S,alpha)
	reltol  = sqrt(eps);
	maxiter = 100;

	p       = 0.7;
	g       = 9.81; % m/s
	rho     = 1000; % kg/m^3
	kappa   = 0.41;
	
	beta = 1;

	% assume regular grid
	n  = length(N);
	dn = N(2)-N(1);
	I  = speye(n);
	D1 = derivative_1_1d(dn,n);

	% right hand side
	b = rho*g*h*S + alpha*N/(N(end)-N(1));
	% initial condition (ignore turbulence term)
	q = sqrt( (beta.*g./(C.^2.*h.^2)) ./ b);
	q([1,end]) = 0;

	% iteratively solve for q
	for idx=1:maxiter
		% eddy viscosity
		% nu_t = kappa/6 u_* h = kappa (sqrt(g) / C) u h 
		nu_t = 1/6*kappa*(sqrt(g)./C)*abs(q);
		nu_t = diag(sparse(nu_t));
%		nu_t = zeros(size(nu_t));

		% set up matrix
		A = diag(sparse((beta.*g./(C.^2.*h.^2)).*q)) - D1*(nu_t*D1);

		% apply boundary conditions
		% the bc for the second element is not koscher, but as function decays to zero near bank this works
		b([1 end]) = 0;
		A(1,:)     = 0;
		A(1,1)     = 1;
%		A(2,1)     = 0;
%		A(end-1,end) = 0;
		A(end,:)   = 0;
		A(end,end) = 1;
%		spy(A)
%		pause

		% solve linear equation
		q_ = A \ b;
		% change
		delta = norm(q_-q)/norm(q_)
%		plot([q q_])
%		pause(0.1)
		q = p*q_ + (1-p)*q;
		% test for convergence
		if (delta < reltol)
			break;	
		end % if delta
	end % for idx
	if (delta > reltol)
		warning('no convergence');
	end % if delta
	% determine velocity
	u = q./h;
end % lateral_division_method

