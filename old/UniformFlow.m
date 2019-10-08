% Wed  3 May 13:55:46 CEST 2017
% returns requested valued based on the definitions
% This would be a perfect example in prolog
% 	Q   = U H W
% 	U   = C sqrt(H S)
%	us  = sqrt(gHS)
%	tau = rho us^2
% TODO, cd, f, manning
classdef UniformFlow < handle
	properties
		% discharge
		Q_
		% flow velocity
		U_
		% shear velocity
		us_
		% shear stress
		tau_
		% channel conductivity (chezy)
		C_
		% channel width
		W_
		% channel depth
		H_
	end

	% TODO, warn if overdetermined
	function obj = UniformFlow(varargin)
		for idx=1:2:length(varargin)
			field = varargin{idx};
			val   = vararing{idx+1};
			if (~isempty(val))
				switch (field)
				case {'C','Q','W','S','U','H','us','tau'}
					obj.([field,'_']) = val;
				otherwise
					error('unknown field');
				end
			
		end
	end
	function U = U(obj)
		U = obj.U_;
		if (isempty(U))
			U = sqrt(obj.H*obj.S);
		
	end

	if (isempty(C))
		if (isempty(S) || isempty(U) || isempty(H))
			error('Require S, U amd H to compute C');
		end
		C = U/sqrt(H*S);
	end
	if (~isempty(Qw) && isempty(U))
		U = (C^2*Q/W*S)^(1/3);
	end
	if (isempty(U))
		if (isempty(S))
			error('either S or U must be provided')
		end
		U = C*sqrt(H*S);
	end
end % classdef UniformFlow

