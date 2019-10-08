% Mon Mar 16 09:59:33 CET 2015
% Karl Kastner, Berlin
% Stationary Rating Curve, Keulegan Roughness formulation
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
classdef KeuleganRatingCurve < RatingCurve
	properties
		opt
		param0
		cext
		cflag
	end % properties
	methods
		function obj = KeuleganRatingCurve(varargin)
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
		function [param obj] = fit_(obj,t0,q0,h0) %,param0)
			% lsqnonlin requires double input
			h0 = double(h0);
			q0 = double(q0);

			if (nargin() < 5)
				param0 = [1e-3; 1e-5];
			end
			[param resn res flag] = lsqnonlin(@(c) obj.objective(c,t0,q0,h0), param0, [0 0], [],obj.opt);
		end % fit_

		function [q g obj] = rcfunc_(obj, c, t, h)
			% (pseudo)-area
			A = double(obj.afunc(h));
			P = double(obj.pfunc(h));
			R = A./P;

			% roughness length
			z0 = c(1,:);
			
			% stationary energy slope
			S0 = c(2,:);

			% depth dependent Chezy's coefficient
			C = sqrt(Constant.g)/Constant.Karman*log(h./(Constant.Euler*z0));
			
			% stationary velocity
			u0 = C.*sqrt(R.*S0);

			% discharge
			q = A.*u0;
			
			% value of the rating curve
%			f = sqrt(Constant.g)/Constant.Karman*log(h./(Constant.Euler*c(1))).*A.*sqrt(c(2).*h);
			
			% gradient of the rating curve value with respect to the parameters
			if (nargout() > 1)
%				g = grad(@(c) obj.rcfunc_(c, t, h), c);
				g = [ sqrt(Constant.g)/Constant.Karman*(-1./c(1)).*A.*sqrt(c(2).*R), ...
				      0.5*sqrt(Constant.g)/Constant.Karman*log(h./(Constant.Euler*c(1))).*A.*sqrt(R./c(2)) ];
			end % nargout
		end % rcfunc_
	end % methods
end % classdef KeuleganRatingCurve

