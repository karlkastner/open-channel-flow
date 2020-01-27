% Sun 19 Jan 13:15:02 +08 2020
function [x, u, rmse] = transverse_velocity_profile_tidal_channel(C, harg, w, S, lambda, n, bc)

% rho g h S + cd u^2 + 1/2 rho lambda sqrt(f/8) d/dn(h^2 d/dn u) = 0
% deep channel, ignore advective acceleration
% d/dt z + duh/dx, h const, dz/dt = -du/dx
% du/dt + g h S + cd u^2 = 0
% S = S_0 + S_1 sin(o t)
% u = u_0 + u_1 sin(o t + phi_1) + u_2 sin(2 o t + phi_2)
% tidally dominant:
% u^2|_0 = u_0^2 + 1/2 \hat u_1^2 + 1/2 \hat u_2^2
% u^2|_1 = 2 u_0 u_1
% u^2|_2 = 1/2 \hat u_1^2 + u_0 u_2
% du/dt = sum i k o u_i

% iteration required, as first powers are required

% TODO tidal dominance

% rho g h S + cd u^2 + 1/2 rho lambda sqrt(f/8) d/dn(h^2 d/dn u) = 0

% compute u^3 profile for tidal or nontidal flow (more concentrated in centre?)
%D2(h^2 u^2) = D2(h^2 (u0^2 + 1/2 u1^2)
	abstol = 1e-7;
	maxiter = 100;
	p = 0.9;

	if (nargin() < 7 || isempty(bc))
		bc = 'dirichlet';
	end
	if (nargin() < 6 || isempty(n))
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
	I  = speye(n);

	o = 2*pi/(86400);

	% momentum for u0
	b = -[(rho*g*h*S(1));
	      (rho*g*h*S(2));
	      (rho*g*h*S(3))];

	% initial condition (no diffusion or interaction
	u = [   0.5+randn(n,3) ] + [zeros(n,1), 1i*randn(n,2)];
	%	zeros(n,1) ];

	Z = sparse(n,n);

	% Note, his is actually rho cd
	cd = -1/8*rho*f;
	l  = 0.5*rho*lambda*sqrt(1/8*f);
% u0 :     cd*u0^2 + l*D2*h^2*u0^2 + cd 1/2 |u1|^2 + l*D2*h^2*(1/2)*|u1|^2
% u1 :                             + cd 2 u0 u1 + l*D2*h^2*2*u0*u1
% u2 :                             + cd 1/2 ..                            m cd 2 u0 u2
% TODO nu missing
	k = 0;
	while (1)
	k = k+1;
	uold = u(:);
A = [ cd*diag(u(:,1)) + l*D2*diag(h.^2.*u(:,1)), cd*1/2*diag(conj(u(:,2))) + l*D2*diag(h.^2.*conj(u(:,2))), cd*1/2*conj(u(:,3)) + l*D2*diag(h.^2.*conj(u(:,3))); ...
                                Z, rho*1i*o*I + 2*cd*diag(u(:,1))         + l*D2*diag(h.^2.*2.*u(:,1)),     Z;
                                Z, 1/2*cd*conj(u(:,2)) + l*D2*diag(h.^2.*conj(u(:,2))), rho*2i*o*I + 2*cd*u(:,1)         + l*D2*diag(h.^2.*2.*u(:,1))];

	switch (bc)
	case {'dirichlet'}
		for jdx=1:3
			A(1+n*(jdx-1),:) = 0;
			A(n*jdx,:)       = 0;
			A(1+n*(jdx-1),1+n*(jdx-1)) = 1;
			A(  n*jdx,n*jdx)           = 1;
			b(1+n*(jdx-1)) = 0;
			b(n*jdx) = 0;
		end
	%case {'neumann'}
	otherwise
		error('here');
	end

	p = 0.25;
	u = p*(A \ b) + (1-p)*uold;
	rmse = rms(u-uold);
	disp(rmse)
%	u0 = u(:,1);
%	u1 = u(:,2);
	u = reshape(u,n,3);
	if (rmse < abstol)
		break;
	end
	if (k>=maxiter)
		warning('no convergence');
		break;
	end
	end % while 1
end % function

