% Sun May  1 17:46:08 MSD 2011
% Karl Kästner
%% set surface level as Dirichlet boundary condition
% TODO, this is actually depth not level, avoid passing of w0
function [A, B, rhs, obj] = bc_level(obj, t, q1, dt, hfun, w0)

		if (obj.QAflag)
			% A
			val = w0*hfun(t);
		else
			% h
			val = hfun(t);
		end
		
		A = [1, 0;	% h(0) = h0
                     0, 0];	%
		B = [0, 0;	% 
                     0, 1];	% dq/dx|0 = 0
		rhs = [val;
		        0];
end % bc_level

