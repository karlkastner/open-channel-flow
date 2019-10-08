% Mo 6. Jul 15:05:15 CEST 2015
% Karl Kastner, Berlin
%
%% stationary 1D shallow water equation across a river section
%% 0 = - g h S0 - tau_b/rho + d/dn (nu h du/dn)
%% 0 = - g h S0 + g u^2/C^2 + d/dn (nu h du/dn)
%
%% includes tranvese gradient term
%% 
%% note that shiono/knight 1991 provide an _analytic_ solution,
%% which takes the form of an expontially decaying side wall effect
%%
% -> why is this not solved in u^2?
function u = stationary_1d_swe(N,h,C,S,alpha,s_nu)
	reltol  = sqrt(eps);
	maxiter = length(N);

	if (nargin()<5)
		s_nu = 1;
	end

	g       = 9.81; % m/s
	rho     = 1000; % kg/m^3
	kappa   = 0.41;
	% relaxation parameter
	p       = 0.65;
	
	beta = 1;

	% assume regular grid
	n  = length(N);
	dn = N(2)-N(1);
	I  = speye(n);
	D1 = derivative_1_1d(dn,n);

	% right hand side
	% energy slope + momentum redistribution by curvature
	b = g*h*S + vander_1d(N/(N(end)-N(1)),1)*alpha;

	% initial condition (ignore turbulence term)
	u = sign(b).*(sqrt( (beta.*g./C.^2) ./ abs(b)));
	% dirichlet bc
	u([1,end]) = 0;
	u(0==b) = 0;

	% iteratively solve for q
	for idx=1:maxiter
		% eddy viscosity
		% nu_t = kappa/6 u_* h = kappa (sqrt(g) / C) u h 
		nu_t = 1/6*kappa*(sqrt(g)./C).*(u.*h);
		nu_t = diag(sparse(nu_t));

		% set up matrix
		A = ( diag(sparse(beta.*g./C.^2.*u)) - D1*(s_nu*nu_t*D1) );

		% apply boundary conditions
		% the bc for the second element is not koscher, but as function decays to zero near bank this works
%		b([1 end]) = 0;
%		A(1,:)     = 0;
%		A(1,1)     = 1;
%		A(end,:)   = 0;
%		A(end,end) = 1;
		switch (bcmode)
		case {'dirichlet'}
		case {'neumann'}
		end

		% solve linear equation
		u_ = A \ b;

		% change
		delta = norm(h.*(u_-u))/norm(h.*u_);

		% step with relaxation
		u = p*u_ + (1-p)*u;

		% test for convergence
		if (delta < reltol)
			break;	
		end % if delta
	end % for idx
	if (delta > reltol)
		warning('no convergence');
		u = NaN(size(u));
	else
		fprintf(1,'converged ind %d steps\n',idx);
	end % if delta
end % function stationary_1d_swe

