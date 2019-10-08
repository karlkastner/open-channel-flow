% Thu 20 Sep 18:28:47 CEST 2018
% Karl Kastner, Berlin
%
%% log profile with linear bend correction
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
classdef Log_profile_with_linear_bend_correction < Log_profile
	properties % (Constant)
		%np = 3;
		%p=1;
	end % properties
	methods (Static)
		df_dc = df_dc_(ln_z0,c,h,s);
		f     = profile_(ln_z0,c,h,s);
	end
	methods
		% constructor
		function obj = Log_profile_with_linear_bend_correction(varargin)
			obj = obj@Log_profile(varargin{:});
			obj.np = 3;
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

		% retrive normalized perturbation parameter
		% pp = (us/kappa)^-1 c = a^-1 c
		function [pp, serr_pp, bias_pp, obj] = perturbation_parameter(obj,pp)
			if (nargin()<2)
			a    = obj.param(1,:);
			c    = obj.param(3,:);
			s2param = obj.s2param;
			sa2  = s2param(1,:);
			sc2  = s2param(3,:);

			sa2(isnan(sa2)) = 0;

			% similar to the roughness length, this parameter is only defined of the
			% shear velocity is known with sufficient accuracy
			fdx      = (sa2 >= a.*a);
			c(fdx)   = NaN;
			sc2(fdx) = NaN;

			% first oder bias corrected perturbation parameter
			% TODO covariance
			pp = c./a - (c./a.^3).*sa2;

			% standard error
			if (nargout() > 1)
	        		serr_pp = sqrt(sc2./(a.*a) + (c./(a.*a)).^2.*sa2);
			end

		        % first order bias of the uncorrected perturbation parameter
			% pp = pp_b - b
			% TODO covariance
			if (nargout() > 2)
			        bias_pp  = (c./a.^3).*sa2;
			end
			else
				a    = obj.param(1,:);
				c = a.*rvec(pp);
				obj.param(3,:)   = c;
				% obj.s2param(3,:) = zeros(size(c));
			end
		end % perturbation_parameter
	end % methods
end % classdef Log_profile_perturbed

