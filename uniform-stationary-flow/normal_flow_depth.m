% 2016-05-17 13:37:38.441944219 +0200
% Karl Kastner, Berlin
%
%% normal flow depth for uniform stationary flow
%% function H = normal_flow_depth(Q,W,C,S)
function H = normal_flow_depth(Q,W,C,S)
	%H = cbrt(Q.^2./(W.^2*S));
	H = (Q.^2./(C.^2.*W.^2.*S)).^(1/3);
end

