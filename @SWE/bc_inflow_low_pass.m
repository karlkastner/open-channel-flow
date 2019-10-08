% Thu Apr 28 04:04:04 MSD 2011
% Karl Kästner
%
%% set low frequency Dirichlet, high frequency pass boundary condition
%	
% h(1) = 0.65; for shock
%
% non-reflecting boundary condition
% extrapolate 0-order
% TODO avoid passing of w0
%
%
function [A B rhs obj] = bc_inflow_low_pass(obj, t, y1, dt, q0fun, Tlp, cd, w0, dzb_dx)
		g = obj.g;

		p = min(1,dt/Tlp);
		q0 = q0fun(t);

		%    q0 = dt/T qin + (1-dt/T)q1, T<dt
		% => q0 = qin + (1-dt/T)dx dq/dx

		if (obj.QAflag)
			A1 = y1(1);
			Q1 = y1(2);
			h1 = A1/w0;
			q1 = Q1/w0;
			q0 = q0/w0;
		else
			h1    = y1(1);
			q1    = y1(2);
		end
		dh_dx = -abs(q1)*q1*cd/g*h1.^-3;

		% TODO, dh/dx = Sf-Sb-Sw (!), not zero
		A = [0,0;	 % dh/dx = 0
                     0,1];	 % q = q0

		B = [1,0;	 % dh/dx = 0
                     0,0];	 % q = q0

		rhs = [w0*p*dh_dx + w0*(1-p)*dh_dx;   % dh/dx = 0
		          w0*p*q0 + w0*(1-p)*q1]; % q = q0
end % bc_inflow_low_pass

