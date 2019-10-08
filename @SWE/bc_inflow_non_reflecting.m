%Mon 20 Nov 15:05:33 CET 2017
%% set non-reflecting boundary condition
function [A,B,rhs,obj] = bc_inflow_non_reflecting(obj, t,q1,dt,q0fun)
        %g = 9.81;                                                               
	% TODO better h0 = h1 + dx(Sf-Sb) for stationary inflow
	q0 = [q1(1);     % dh/dx = 0
              q0fun(t)]; % q = q0

	[L R Ri] = obj.fluxmateig(q1);

	% leaving wave
	m = diag(L)<0;

	% entering wave
	p = diag(L)>0;

	rhs =   R*diag(m)*Ri*q1 ...
	      + R*diag(p)*Ri*q0;

	A = [1 0;
	     0 1];

	B = [0 0;
	     0 0];	
end


