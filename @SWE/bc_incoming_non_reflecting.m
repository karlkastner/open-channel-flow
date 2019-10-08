%Mon 20 Nov 15:05:33 CET 2017
% TODO avoid passing of w0, pass id1
%% set non-reflecting boundary condition for the 1D SWE
function [A,B,rhs,obj] = bc_incoming_non_reflecting(obj, t, q1, dt, q0fun, w0)
	q0 = q0fun(t);

	[L R Ri] = swe.fluxmateig(q1,w0);

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


