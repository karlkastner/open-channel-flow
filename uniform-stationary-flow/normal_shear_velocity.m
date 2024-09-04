% Wed  1 Jan 16:39:36 +08 2020
% function us = normal_shear_velocity(Q,W,C,S,ismanning)
function us = normal_shear_velocity(Q,W,C,S,ismanning)
	if (nargin<5)
		ismanning = false;
	end
	if (issym(Q))
		syms g
	else
	g = Constant.gravity;
	end
	H = normal_flow_depth(Q,W,C,S,ismanning);
	U = Q./(H*W);
	if (ismanning)
		C = manning2chezy(C,H);
	end
	us = sqrt(g)./C.*U;
end
