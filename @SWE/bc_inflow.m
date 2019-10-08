% Thu Apr 28 04:04:04 MSD 2011
% Karl Kästner
	
%% inflow boundary condition
% TODO accept variable width and depth at last interface dh/dx = Sf-S0
function [A, B, rhs, obj] = bc_inflow(obj, t, y1, dt, Q0fun, cd, w0, dzb_dx)
		g = obj.g;

		A = [0,0;	 % dh/dx = 0
                     0,1];	 % q = q0

		B = [1,0;	 % dh/dx = Sf
                     0,0];	 % q = q0

		if (obj.QAflag)
			A1 = y1(1);
			Q1 = y1(2);	
			h1 = A1/w0;
			q1 = Q1/w0;
			Q0 = Q0fun(t);
		else
			h1 = y1(1);
			q1 = y1(2);
			Q0 = Q0fun(t)/w0;
		end

		% TODO maybe it would best to simply linearly extrapolate the area
		q0    = Q0/w0;
		q     = q1; %0.5*(q0+q1);	
		dh_dx = -abs(q)*q*cd/g*h1.^-3;
		%dh_dx = -abs(q1)*q1*cd/g*h1.^-3;

		if (obj.QAflag)
			rhs   = [w0*dh_dx; % dA/dx = w*Sf
			         Q0]; % q = q0
		else
			rhs   = [   dh_dx; % dh/dx = Sf
			            q0]; % q = q0
		end
end % bc_inflow

