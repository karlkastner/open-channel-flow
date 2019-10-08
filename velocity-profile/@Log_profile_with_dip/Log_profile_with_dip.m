% Sat  7 Jan 15:52:46 CET 2017
% Karl Kastner, Berlin
%
%% Logarithmic profile with dip
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
classdef Log_profile_with_dip < Log_profile
	properties (Constant)
		np = 3;
	end % properties
	methods
		% constructor
		function obj = Log_profile_with_dip(param,sparam)
			if (nargin() > 0)
				obj.param = param;
			end
			if (nargin() > 2)
				obj.sparam = sparam;
			else
				obj.sparam = NaN(size(param));
			end
		end % Log_profile

		% in external file
		% function fit(U,S,mask,h)
		% end
		% function predict(S,h)
		% end
		
		% retrieves shear velocity from linear parameters
		% inherited from Log_profile
		% function [us serr_us obj] = shear_velocity(obj)

		% retrives roughness length from linear parameters
		% inherited from log_profile

		% retrive noramlized dip parameter
		% dp = (us/kappa)^-1 c = a^-1 c
		function [dp, serr_dp, bias_dp, obj] = dip_parameter(obj)
			a    = obj.param(1,:);
			c    = obj.param(3,:);
			sa2  = obj.s2param(1,:);
			sc2  = obj.s2param(3,:);

			% as for the log profile, the roughness and the dip parameter
			% can only be estimated if the standard error of the shear velocity
			% is sufficiently low
			fdx      = (sa2 >= a.*a);
			c(fdx)   = NaN;
			sc2(fdx) = NaN;

			% dip parameter
			% TODO first order bias correction
			dp = c./a;

			% standard error of dip parameter
			% ln_z0.serr = sqrt(sb2./(a.*a) + (b./(a.*a)).^2.*sa2 + 2*b./a.^3.*sab);
			if (nargout() > 1)
				serr_dp = sqrt(sc2./(a.*a) + (c./(a.*a)).^2.*sa2 + 2*c./a.^3.*sac);
			end
			
			% first order bias of dip parameter
			% dp = dp_b - b
			% TODO
			if (nargout() > 2)
				bias_dp = NaN;
			end
		end % dip_parameter
	end % methods
end % classdef Log_profile_with_dip

