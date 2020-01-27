% 2016-05-17 13:37:38.441944219 +0200
% Karl Kastner, Berlin
%
%% normal flow discharge for uniform stationary flow
% function Q = normal_flow_discharge(H,W,C,S)
function Q = normal_flow_discharge(H,W,C,S,ismanning)
	if (nargin()>4 && ismanning)
		% C arg is actually n
		C = manning2chezy(C,H);
	end
	%H = cbrt(Q.^2./(W.^2*S));
%	H = (Q.^2./(C.^2.*W.^2.*S)).^(1/3);
	U = C.*sqrt(H.*S);
	A = H.*W;
	Q = A.*U;
end

