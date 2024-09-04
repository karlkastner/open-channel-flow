% Wed 29 Nov 15:52:11 CET 2017
% Tue 20 Mar 13:56:47 CET 2018
% Karl Kastner, Berlin
%% energy slope (surface slope) for uniform stationary flow
%
%% normal flow slope in uniform stationary flow
% function S = normal_flow_slope(Q,H,W,C)
function S = normal_flow_slope(Q,H,W,cf,type)
	if (nargin()<5 | isempty(type))
		type = 'chezy';
	end
	switch (lower(type))
	case {'drag','cd'}
		Cz = drag2chezy(cf);
	case {'c','cz','chezy'}
		Cz = cf;
	case {'manning'}
		Cz = manning2chezy(cf,H);
		%cf.*H.^(1/6);
	otherwise
		error('here');
	end
	S = sign(Q).*1./H.^3.*(Q./(Cz.*W)).^2;
end

