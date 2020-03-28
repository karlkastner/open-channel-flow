% Thu Apr 28 04:50:48 MSD 2011
% Karl Kästner
%
%% set non-reflecting boundary condition
%% extrapolate 0-order
function [A, B, rhs, obj] = bc_nonreflecting(obj,t,q,dt)
		A = [ 0 0;
		      0 0];
		B = [ 1,0;	% dh/dx = 0
		      0,1];	% dq/dx =0
		rhs = [0;0];
end % bc_nonreflectig

