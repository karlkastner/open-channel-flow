% Sun 23 May 10:27:07 CEST 2021
% Karl Kastner, Berlin
%
%% roughness coefficient from uniform stationary flow
%function rgh = normal_flow_roughness(Q,H,W,S,type)
function rgh = normal_flow_roughness(Q,H,W,S,type)
	U = Q./(H.*W);
	rgh = U./sqrt(H.*S);	

	if (nargin()>4)
	switch (lower(type))
	case {'n','manning'}
		rgh = chezy2drag(C,H);
	case {'chezy','cz'}
		% nothing to do
	case {'drag','cd'}
		rgh = chezy2drag(C);
	otherwise
		error('here');
	end
	end
end

