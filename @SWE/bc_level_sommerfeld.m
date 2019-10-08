% Mon 20 Nov 12:21:19 CET 2017
% Karl Kästner
%% set surface level as boundary condition by sommerfeld method
function [A B rhs obj] = bc_level_sommerfeld(obj, t, q, dt, hfun)
		g = obj,g;
		h0 = hfun(t);
		% TODO, how to choose c?
		c = sqrt(g*h0);
		A = [1, 0;	% h(0) = h0
                     0, 0];	% sommerfeld
		B = [0, 0;	% h0 = 0
		     sqrt(g/h0), -1];
                %     c,-1];	% sommerfeld: c dh/dx - dq/dx = 0
		%     0 = qr-ql + 2 sqrt(gh)(sqrt hl - sqrt(hr))
		rhs = [h0;
		       0];
end % bc_level

