% 2016-05-17 13:37:38.441944219 +0200
% Karl Kastner, Berlin
%
%% normal flow depth for uniform stationary flow
%% function H = normal_flow_depth(Q,W,C,S)
function H = normal_flow_depth(Q,W,cf,S,ismanning)
	%H = cbrt(Q.^2./(W.^2*S));
	if (nargin()>4 && ~isempty(ismanning) && ismanning)
		n = cf;
		H = ( (Q.*n)./(sqrt(S).*W) ).^(3/5);
	else
		C = cf;
		H = (Q.^2./(C.^2.*W.^2.*S)).^(1/3);
	end
end

