% 2017-01-09 19:14:31.475367378 +0100
% Karl Kastner, Berlin
%
%% propagation of a diffusion wave (flood wave), c.f. ponce
%
function [tt, x, QQ] = diffusion_wave(T,L,Q0func,bcfunc,C,S0,W,dt,dx,dfactor,method,tfactor)

	if (nargin()<12)
		tfactor = 0;
	end

%method = 'crank-nicolson';
%method = 'upwind';

% discretise in space
nx  = round(L/dx);
dx  = L/(nx-1);
x   = dx*(0:nx-1)';

switch (method)
case {'crank-nicolson','btcs'}
	D1x = derivative_matrix_1_1d(nx,L);
case {'upwind'}
	if (dt > 0)
		D1x = derivative_matrix_1_1d(nx,L,1);
	else
		D1x = derivative_matrix_1_1d(nx,L,-1);
	end
otherwise
	error('here');
end
%D1xr = derivative_matrix_1_1d(nx,L,1);
%D1xr(1,1:2) = 1/dx*[-1 1];
%D1xr(end,end-1:end)
%pause
%D1xr(1,
%D1xr(1,:) = 0;
%D1xr(end,:) = 0;

D2x = derivative_matrix_2_1d(nx,L);

% initial condition
Q   = Q0func(x);

% solve
t  = T(1);
QQ = [];
tt = [];
while (true)
	% write back
	QQ(:,end+1) = Q;
	tt(end+1)   = t;

	% TODO check cfl
	if (dt > 0)
		% forward in time
		% limit last time step
		dt = min(dt,T(2)-t);
		if (dt <= 0)
			break;
		end
	else
		dt = max(dt,T(2)-t);
		if (dt>=0)
			break;
		end
	end

	% step
	[Q t]= step(t,Q,C,W,S0);

	% store
%	QQ(:,end+1) = Q;
%	tt(end+1)   = t;
end % while true

function [Q t] = step(t,Q,C,W,S0)
	%Ponce
	%dQ/dt + cdQ/dx - d d2^Q/dx2 = 0
	%d = Q/(2WS0)
	%c = dQ/dA = 3/2u
	
	% flow depth (chezy, statinary flow in wide rectangular channel)
	h0 = (Q/(W*C*sqrt(abs(S0)))).^(2/3);
	% area
	A0 = W*h0;
	% velocity
	u0 = Q./A0;
	% celerity
	c0 = 3/2*u0;
	% diffusion coefficient
	d = 0.5*Q/(W*abs(S0));
	d = dfactor*d;

	% cfl    = max(abs(c0))*abs(dt/dx);
%	dt_max = 0.8*(dx/max(abs(c0)));

	
	% time derivative matrix
	A1 = -diag(sparse(c0))*D1x;
	A2 =  diag(sparse(d))*D2x;

	% step in time
	Q__ = Q;

	% apply bc
	% TODO apply bc by bcfunc
	Q_bc = bcfunc(t+dt,Q);
	I  = speye(size(A1));
	if (dt > 0)
		% TODO this seems not quite right, because the new value is then used to compute the derivative at the second row
		Q(1)    = Q_bc(1);
		A1(1,:) = 0;
		A2(1,:) = 0;

		% TODO extrapolate coefficients as well
		if (~strcmp(method,'upwind'))
			% for upwind, the value n+1 is not used and does not need to be extrapolated
			A1(end,end-1:end) = A1(end,end-1:end) + (-0.5*c0(end)/dx)*[-1 2]; % makes derivative first order
		end
		A2(end,end-1:end) = A2(end,end-1:end) + (d(end)/dx^2)*[-1 2]; % effectively sets diffusion to zero
		%A(end,end-1:end) = A(end,end-1:end) + (-0.5*c0(end)/dx + d(end)/dx^2)*[-1 2];
	else
		Q(end)    = Q_bc(end);
		A1(end,:) = 0;
		A2(end,:) = 0;

		if (~strcmp(method,'upwind'))
			% for upwind, the value -1 is not used and does not need to be extrapolated
			A1(1,1:2) = A1(1,1:2) + (+0.5*c0(1)/dx)*[2 -1];
		end
		A2(1,1:2) = A2(1,1:2) + (d(1)/dx^2)*[2 -1];

%		I(1,1:3)  = [0 2 -1];
%		A1(1,1:4) = 2*A1(2,1:4) + -1*A1(3,1:4);
%		A2(1,1:4) = 2*A2(2,1:4) + -1*A2(3,1:4);
%		A1(1,1:2) = A1(1,1:2) + (+0.5*c0(1)/dx)*[2 -1];
%		A2(1,1:2) = A2(1,1:2) + (d(1)/dx^2)*[2 -1];
		% A(1,1:2) = A(1,1:2) + (+0.5*c0(1)/dx + d(1)/dx^2)*[2 -1];
	end

	% step in time
	switch (method)
	case {'upwind'}
		% limit time step
		dt_max = 0.8*min(1./(abs(c0)/dx + 2*abs(d)/dx^2));
		dt = sign(dt)*min(abs(dt),dt_max);

		if (0==tfactor)
			A = A1+A2;
			Q = (I + dt*A)*Q;
			%% advection
			%Q = (I+dt*A1)*Q;
			%% diffusion
			%Q = (I+dt*A2)*Q;
		else
			% operator splitting
			% advection
			Q = (I+dt*A1)*Q;
			% diffusion
			% Q = (I+dt*A2)*Q;
			
			%A = (I+dt*A2);
			%Q = Q + (A'*A + tfactor*D'*D) \ (A'*Q);
			dQ = dt*A2*Q;
			% regularise
			D = D2x;
			D(1,1:3) = D(2,1:3);
			D(end,end-2:end) = D(end-1,end-2:end);
			dQ = (I + tfactor*dx^4*D'*D) \ dQ;
%			eigs(inv((I + tfactor*dx^2*D'*D)),4,1e-7)
%pause
%			norm(dQ)
			Q = Q + dQ;

			%Q = Q + dt*((A2'*A2 + tfactor*D'*D)\(A'*Q));
			%A = (I + dt*A2);
			%Q = (A'*A + tfactor*D2x'*D2x) \ (A'*Q);
			% = (A2'-(A2'*A2 + tfactor*D2'*D2)^-1)\(A2'*Q)
		end
		%Q = bcfunc(t+dt,Q);
		% apply boundary condition
		% Q = bcfunc(t+dt,Q);
	case {'btcs'}
		% unconditionally stable
		if (0 == tfactor)
			A = A1+A2;
			Q = (I-dt*A) \ Q;
		else
			% operator splitting
			% advection
			Q = (I-dt*A1) \ Q;
			% diffusion
			if (1)
				D = dx*D1x;
				%D = dx^2*D2x;
				D(1,:) = 0;
				D(end,:) = 0;
				%D(1,1:3) = D(2,1:3);
				%D(end,end-2:end) = D(end-1,end-2:end);
				%dQ = (I-dt*A2)\Q - Q;
				%dQ = (I + tfactor*dx^4*D'*D) \ dQ;
				%Q  = Q + dQ;
				A = (I-dt*A2);
				Q = (A'*A + tfactor*D'*D) \ (A'*Q);
			else
				A = (I-dt*A2);
				[u s v] = svd(full(A));
				s = max(0.5,diag(s));
				Q = (v'*diag(sparse(1./s))*u')*Q;
			end
			%Q = (I-dt*A2) \ Q;
			%dQ = (I - (I-dt*A2)^-1) \ Q
			%Q = (A2'-(A2'*A2 + tfactor*D2'*D2)^-1)\(A2'*Q)
			% A = A1+(A2'*A2 + tfactor*D2'*D2)\A2'
		end
		
	case {'crank-nicolson'}
		if (0 == tfactor)
			A = A1+A2;
			Q = (I-0.5*dt*A) \ (Q+(0.5*dt)*(A*Q));
		else
			% operator splitting : advection
			Q = (I-dt*A1) \ Q;
			
		end
	otherwise
		error('here');
	end
%		% thikonov regularisation (for inverse problems)
%		if (dt < 0 && tfactor > 0)
%			dQ = Q-Q__;
%			dQ_ = dQ;
%			dQ = (I+tfactor*dt*dt*(D1xr'*D1xr)) \ dQ;
%			dQ = dQ+(mean(dQ_)-mean(dQ));
%			Q  = Q__+dQ; 
%if (0)
%			figure(2)
%			clf
%			subplot(2,2,1)
%			plot([dQ dQ_ dt*A*Q__])
%			subplot(2,2,2)
%			plot([Q__ Q])
%			pause(0.1)
%end
%		end

	%Q = Q + dt*dQ_dt;
	t = t+dt;

	% apply boundary condition afterwards
end % dQ_dt;


%beta 0.01
%Se = Ke/2g d(Q/A)^2/dx
%% where is the bed slope?
%Wf : wind shear
%q : inflow
%B : top width
%
%% friction slope
%Sf = 1/C^2 V^2/R;
%
%% eddy slope
%Se = 0;
%
%% chow 1988
%% d(A+A0)/dt + dQ/dx = q
%% dQ/dt + d/dx betaQ^2/A + gA(dh/dx + Sf + Se) - beta q_i v_i + Wf B = 0
%
%% A0 ignored
%Adot = Dx*Q;
%
%% inflow and wind shear ignored
%Qdot = Dx*(beta*(Q.^2./A) + g*A) + g*A*(Sf + Se)
%
%dot = [Adot; Qdot];
%
% dQ/dx + 1/c dQ/dQ = Qin
% dQ/dx + alpha beta Q^(beta-1) dQ/dt = q

end % diffusion wave

