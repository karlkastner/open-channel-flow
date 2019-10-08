% Mon Mar 16 11:33:29 CET 2015
% Karl Kastner, Berlin
%
%% Dynamic Rating Curve, Keulegan roughness formulation
%% (dynamic = correction for hysteresis loop)
% implmentation assumes for the time being that the hyrdaulic radius is H
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
classdef DynamicKeuleganRC < RatingCurve
	properties
		opt
		param0
		cext
		cflag
	end % properties
	methods
		function obj = DynamicKeuleganRC(varargin)
			obj.cflag = false(3,1);
			for idx=1:2:length(varargin)
				switch(varargin{idx})
				case {'afunc'}
					obj.afunc = varargin{idx+1};
				case {'pfunc'}
					obj.pfunc = varargin{idx+1};
				otherwise
					error('PowerRatingCurve','Unknown Parameter');
				end % switch
			end % for
			obj.rcfunc = @obj.rcfunc_;
			obj.opt = optimset( 'MaxIter', 10000, ...
				'MaxFunEvals', 10000, ...
				'Algorithm', 'levenberg-marquardt', ...
				'Jacobian', 'on', ...
				'TolFun', 1e-12, ...
				'TolX', 1e-12, ...
				'TolCon', 1e-12);
		end % PowerRatingCurve
		function [param obj] = fit_(obj,t0,q0,h0, dh_dt0)
			% lsqnonlin requires double input
			h0 = double(h0);
			q0 = double(q0);
			dh_dt0 = double(dh_dt0);
			if (nargin() < 6)
				param0 = [1e-3; 1e-5];
			end
			[param resn res flag] = lsqnonlin(@(c) obj.objective(c,t0,q0,h0,dh_dt0), param0, [0 0], [],obj.opt);
		end % fit_

		function [q g obj] = rcfunc_(obj, c, t, h, dh_dt)
			% (pseudo)-area
			A = double(obj.afunc(h));
			P = double(obj.pfunc(h));
			R = A./P;

			% roughness length
			z0 = c(1,:);

			% stationary energy slope
			S0 = c(2,:);

			% depth dependent Chezy coefficient
			C = sqrt(Constant.g)/Constant.Karman*log(h./(Constant.Euler*z0));

			% stationary velocity
			u0 = C.*sqrt(R.*S0);

			% instationary velocity
			u  = u0.*sqrt(1 + min(1,max(-1, 2./(3*u0*S0).*dh_dt)));

			% instationary discharge
			q  = A.*u;

			% gradient of the rating curve value with respect to the parameters
			if (nargout() > 1)
				g = grad(@(c) obj.rcfunc_(c, t, h, dh_dt), c);
%				g = [ sqrt(Constant.g)/Constant.Karman*(-1./c(1)).*A.*sqrt(c(2).*h), ...
%				      0.5*sqrt(Constant.g)/Constant.Karman*log(h./(Constant.Euler*c(1))).*A.*sqrt(h./c(2)) ];
			end % nargout
		end % rcfunc_
	end % methods
end % classdef DynamicKeuleganRC

