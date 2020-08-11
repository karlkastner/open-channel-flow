% Thu 16 Apr 16:10:23 +08 2020
% Thu 23 Apr 15:36:31 +08 2020
classdef Cross_Section < handle
	properties
		% nt x nn
		n = struct( ...
			... % across-channel coordinates of integration points
			'x', [], ...
			... % specific discharge
			'q', [], ...
			... % chezy
			'C', [], ...
			... % median grain size
			'd50', [], ...
			... % 90th-percentile
			'd90', [], ...
			... %dsd
			... % bed level
			'zb', [], ...
			... % integration weights
			'weight', [] ...
			);
		% integrated quantities
		%i = struct ( ...
		% [ntx1] surface elevation
		zs = [];
		zb = struct();

		% [ntx1] area
		discharge
		area  = [];
		Qs

		% [ntx1] surface slope
		slope  = [];
	
		% TODO mae dxx, pxx	
		d50
		d90
		c

		% [scalar] top width
		width
		width_

		% temperature
		T_C

		p

		nsmooth

		method
		nn
		order
		flag
		verbose
		zblim
		zslim
		dparam_
		mode
		widthadapt
		force_ordered = false;
		dsd
	end
	methods
		function obj = Cross_Section(varargin)
			obj = parse_arguments(obj, varargin{:});
		end

		function dparam = dparam(obj)
			if (isempty(obj.dparam_))
				obj.dparam_ = lognfit_quantile([0.5,0.9],[obj.d50,obj.d90]);
			end
			dparam = obj.dparam_;
		end
	end
end % Cross_Section

