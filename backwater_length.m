% 2022-01-04 19:51:14.528064241 +0100
function L = backwater_length(Q,W,h0,C,Sb,varargin)
	hinf = normal_flow_depth(Q,W,C,Sb,varargin{:})
	S0   = normal_flow_slope(Q,h0,W,C,varargin{:})
	L = -(h0-hinf)./(S0-Sb);
end

