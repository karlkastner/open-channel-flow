% Wed 29 Nov 15:52:11 CET 2017
% Karl Kastner, Berlin
%
%% friction slope (surface slope) for uniform stationary flow
% function S = friction_slope(Q,H,W,C)
function S = friction_slope(Q,H,W,C)
	%H = cbrt(Q.^2./(W.^2*S));
%	H = (Q.^2./(C.^2.*W.^2.*S)).^(1/3);
	
	% Q^2 = W C H^3 S
	S = Q.^2./(W.^2.*C.^2.*H.^3);
end

