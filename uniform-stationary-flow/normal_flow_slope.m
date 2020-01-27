% Tue 20 Mar 13:56:47 CET 2018
% Karl Kastner, Berlin
%
%% normal flow slope in uniform stationary flow
% function S = normal_flow_slope(Q,W,H,C)
function S = normal_flow_slope(Q,W,H,C,ismanning)
	if (nargin()>4&ismanning)
		C = C*H^(1/6);
	end
	S = 1./H.*(Q./(H.*C.*W)).^2;
end

