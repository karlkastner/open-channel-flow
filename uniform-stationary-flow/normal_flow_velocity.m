% 2017-09-14 13:37:08.868120508 +0200
%
%% normal flow velocity in uniform stationary flow
% function U = normal_flow_velocity(Q,W,C,S)
%
function U = normal_flow_velocity(Q,W,C,S,ismanning)
	if (nargin<5)
		ismanning = [];
	end
	H = normal_flow_depth(Q,W,C,S,ismanning);
	U = Q./(H*W);
%	if (nargin()>4)
%		n = C;
%		C = manning2chezy(n,H);
%	end
%	U = C.*sqrt(H.*S);
end

