% Sat  7 Jan 14:52:28 CET 2017
%
%% logarithmic profile of the streamwise velocity
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
classdef Log_profile < Vertical_profile
	properties %(SetAcces = protected) %(Constant)
		np = 2;
	end % properties
	methods (Static)
		df = profile_(ln_z0,h,s);
		%df = df_dln_z0_(H,z,ln_z0);
		%df = df_dh_(H,hi,ln_z0);
		df = df_dln_z0_(ln_z0,h,s);
		df = df_dh_(ln_z0,h,s);
	end
	methods
		% constructor
		function obj = Log_profile(varargin)
			obj = obj@Vertical_profile(varargin{:});
		end % Log_profile

		% in external file
		% function fit(U,S,mask,h)
		% end
		% function predict(S,h)
		% end
		
		% retrieves shear velocity from linear parameters
		function [us, serr_us, obj] = shear_velocity(obj,us)
			if (nargin() < 2) % get
				a    = obj.param(1,:);
	
			        % estimate of u_s (always unbiased)
				%if (issym(a) || issym(sa2))
				%	syms kappa
				%else
				%	kappa = Constant.Karman;
				%end
				us    = obj.karman*a;
				
	        		% standard error
				if (nargout() > 1)
					s2param    = obj.s2param;
					sa2        = s2param(1,:);
				        serr_us    = obj.karman*sqrt(sa2);
				end
			else % set
				%if (issym(us))
				%	obj.param = sym(obj.param);
				%	obj.s2param = sym(obj.s2param);
				%	syms kappa
				%else
				%	kappa = Constant.Karman;
				%end
				obj.param(1,:)   = us/obj.kappa;
				% obj.s2param(1,:) = zeros(1,length(us));
			end
		end % shear_velocity

		% retrives roughness length from linear parameters
		function [z0 obj] = roughness_length(obj,z0)
			if (nargin()>1)
				obj.ln_z0(log(z0));
			else
				z0 = exp(obj.ln_z0());
			end
		end % roughness_length

		% Log of Roughness length: ln_z0 = -b./a + mu_ln_Z + (b./a.^3).*sa2;
		function [ln_z0, serr_ln_z0, bias, obj] = ln_z0(obj,ln_z0)

			if (nargin() < 2) % get
				a       = obj.param(1,:);
				b       = obj.param(2,:);
			
				% first order bias corrected ln_z0 estimate
			        ln_z0 = -b./a;
				
				s2param = obj.s2param;
				sa2     = s2param(1,:);
				sb2     = s2param(2,:);
	
				sa2(isnan(sa2)) = 0;
		
					% convergence for roughness length requires serr_a^2 < a^2,
					% i.e. the absolute value of the shear velocity	must be larger than its standard error
					% which also shows that averaging before reparametrisation is
					% advantageous, as long as the expected value of the
					% shear velocity (a) does not change between samples
					% this does not effect the convergence of the shear velocity,
					% i.e. the shear velocity can always be estimated, whilst the
					% the roughness length cannot
	
				if (~issym(a) && ~issym(sa2))
					fdx      = (sa2 >= a.*a);
					b(fdx)   = NaN;
					sb2(fdx) = NaN;
				end
	
			        % first order bias of the uncorrected estimate
				% lz = lz_b - b
				bias  = -(b./a.^3).*sa2;
				% bias correction
				ln_z0 = ln_z0 - bias;

				% standard error
				if (nargout() > 1)
					% TODO cov is missing
		        		serr_ln_z0 = sqrt(sb2./(a.*a) + (b./(a.*a)).^2.*sa2);
				end
			else % set
				a = obj.param(1,:);
				b = -a.*rvec(ln_z0);
				obj.param(2,:) = b;
				%obj.s2param(2,:) = zeros(1,length(b));
			end
		end % roughness_length
	end % methods
end % Log_profile

