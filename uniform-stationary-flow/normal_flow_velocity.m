% 2017-09-14 13:37:08.868120508 +0200
%
%% normal flow velocity in uniform stationary flow
% function U = normal_flow_velocity(Q,W,C,S)
%
function U = normal_flow_velocity(Q,W,C,S,varargin)
	H = normal_flow_depth(Q,W,C,S,varargin{:});
	U = Q./(H.*W);
end

