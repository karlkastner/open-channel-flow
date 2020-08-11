% 2016-11-25 09:15:25.983167158 +0100
%% transverse profile of the streamwise velocity
%% c.f. shiono knight
% momentum equation
% TODO ignores effect of transversal slope on bed shear stress
%
% du/dt = 0, du/dx = 0
%
% rho g h S + cd u^2 + 1/2 rho lambda sqrt(f/8) d/dn(h^2 d/dn u) = 0
%
function [x, u, u2] = transverse_velocity_profile(C,harg,w,S,lambda,Q,n,bc)
	if (nargin() < 8 || isempty(bc))
		bc = 'dirichlet';
	end
	if (nargin() < 7 || isempty(n))
		n = 100;
	end
	if (nargin() < 5 || isempty(lambda))
		% shono-knight 1991: lambda = 0.07 in channel centre, and larger near the bank
		% higher with strong secondary flow
		% 1989: lambda = 0.45 in centre
		% 1991 : lambda_floodplain ~ l_channel (2 D_r)^-4, D_r = (h_lf-hmc) / h_mc
		% lambda_a : apparent ~ 10 lambda turb-theory
		lambda = 0.5;
	end
	if (nargin()<6)
		Q = [];
	end

	x = w*((0:n-1)'/(n-1)-0.5);
	if (isa(harg,'function_handle'))
		h = harg(x);
	else
		if (isscalar(harg))
			h = harg*ones(n,1);
		else
			h = harg;
		end
	end

	rho = Constant.density.water;
	g   = Constant.gravity;

	%f = sqrt(8*g)/C^2;
	% c.f. Henderson 1966, 4-8
	f = 8*g/C^2;

	% second derivative matrix
	D1 = derivative_matrix_1_1d(n,w);
	D2 = derivative_matrix_2_1d(n,w);
	I  = eye(n);
	
	% turbulent diffusion
	% Ad = 0.5*rho*lambda*diag(sparse(h.^2.*sqrt(1/8*f)))*D2;
	% d/dx(h^2 d/dx u^2) = h^2 d^2/dx^2 u^2 + d/dx h^2 du^2/dx
	Ad = 0.5*rho*lambda*sqrt(1/8*f)*(  diag(sparse(h.^2))*D2 ...
					 + diag(sparse(D1*h.^2))*D1);

	% friction
	Af = - 1/8*rho*f*I;

	A  = Ad + Af;
	
	% pressure, constant term (rhs)
	b  = -(rho*g*h*S).*ones(n,1);

	switch (bc)
	case {'dirichlet'}
		A(1,:) = 0;
		A(end,:) = 0;
		A(1,1) = 1;
		A(end,end) = 1;
		b(1) = 0;
		b(end) = 0;
	case {'neumann'}
		A(1,:) = 0;
		A(1,1:2) = [-1 1];
		A(end,:) = 0;
		A(end,end-1:end) = [-1 1];
		b(1) = 0;
		b(end) = 0;
	otherwise
		error('here');
	end
	
	% solve for U^2
	u2 = A \ b;
	u  = sqrt(u2);

	% scale to full discharge
	% Q = mean(normal_flow_discharge(h,w,C,S));
	if (~isempty(Q))
		u = u*Q/(w*sum(mid(h.*u))/(n-1));
	end

if (0)
	% inner region
	mu = mean(u);
	fdx = u>mu;
%	fdx = (x>-w/2+2*h & x<w/2-2*h);
	Adu2 = Ad*u2;
	Afu2 = Af*u2;
	[norm(Adu2)./norm(b) norm(Afu2)./norm(b)]
%./norm(b) norm(Ad(:,fdx)*u2(fdx))./norm(b(fdx))]
	[norm(Adu2(fdx))./norm(b(fdx)) norm(Afu2(fdx))./norm(b(fdx))]
	w
	clf
	plot(x,u)
	hold on
	plot(x(fdx),u(fdx))
	pause
end

end % transverse_velocity_profile

