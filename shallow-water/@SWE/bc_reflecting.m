% Thu Apr 28 04:03:13 MSD 2011
% Karl Kästner

%% set reflecting boundary condition
%% extrapolate 0-order and invert v
function [A, B, rhs, obj] = bc_reflecting(obj,t,q,dt)
		A = [0,0;	% dh/dx|0 =
		     0,1];	% q|0 = 0
		B = [1,0;	% dh/dx|0 = 0
		     0,0];	% q|0 = 0
		rhs = [0;0];
end % bc_reflecting

