% 2016-05-17 13:37:38.441944219 +0200
% Karl Kastner, Berlin
%
%% normal flow discharge for uniform stationary flow
% function Q = normal_flow_discharge(H,W,C,S)
% TODO, the sign is wrong (!)
function Q = normal_flow_discharge(H,W,C,S,type)
	if (nargin()>4)
	switch (lower(type))
	case {'n','manning'}
		% C arg is actually n
		C = manning2chezy(C,H);
	case {'chezy','cz'}
		% nothing to do, default
	case {'drag','cd'}
		C = drag2chezy(C);
	otherwise
		error('here');
	end
	end
	%H = cbrt(Q.^2./(W.^2*S));
%	H = (Q.^2./(C.^2.*W.^2.*S)).^(1/3);
	U = C.*sqrt(H.*S);
	A = H.*W;
	Q = A.*U;
end

