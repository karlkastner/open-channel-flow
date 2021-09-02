% Mon 24 May 12:45:10 CEST 2021
% Karl Kastner, Berlin
%
%% normal flow width for uniform stationary flow
% function W = normal_flow_discharge(Q,H,C,S,type)
function W = normal_flow_width(Q,H,C,S,type)
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
	A = Q./U;
	W = A./H;
end

