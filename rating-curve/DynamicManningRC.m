% Mon Mar 16 09:59:33 CET 2015
% Karl Kastner, Berlin
%% Dynamic Rating Curve, Manning roughness formulation
%% (dynamic = correction for hysteresis loop)
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
classdef DynamicManningRC < RatingCurve
	properties
		opt
		param0
		cext
		cflag
	end % properties
	methods
		function obj = DynamicManningRC(varargin)
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
				'Jacobian',  'on', ...
				'TolFun',    1e-12, ...
				'TolX',      1e-12, ...
				'TolCon',    1e-12);
		end % constructor
		function [param obj] = fit_(obj,t0,q0,h0, dh_dt0)
			area       = obj.afunc(h0);
			perimeter  = obj.pfunc(h0);
			radius     = area./perimeter;
			param0     = [1e-1; 1e-4];
			% for numerical reasons the first parameter is defined as the inverse of n
			[param resn res flag] = lsqnonlin(@(c) obj.objective(c,t0,q0,h0,dh_dt0), param0, [0 0], [], obj.opt)
		end % fit_

		function [q g obj] = rcfunc_(obj, c, t, h, dh_dt)
			area       = double(obj.afunc(h));
			perimeter  = double(obj.pfunc(h));
			radius     = area./perimeter;

			% Chezy's coefficient
			C          = c(1,:)*radius.^(1/6);

			% stationary velocity
			u0         = C.*sqrt(radius*c(2,:));
			
			% discharge
			q          = area.*u0.*sqrt(1 + 1./(1.5*u0*c(2,:)).*dh_dt);

			% gradient of the rating curve value with respect to the parameters
			if (nargout() > 1)
				g  = grad(@(c) obj.rcfunc_(c, t, h, dh_dt), c);
%				g  = [ -c(1)^2*A.*h.^(2/3)*sqrt(abs(c(2))).*sign(c(2)), ...
 %                                     0.5*c(1)*A.*h.^(2/3)*1./sqrt(abs(c(2))).*sign(c(2))];
			end % nargout
		end % rcfunc_
	end % methods
end % classdef DynamicManningRC

