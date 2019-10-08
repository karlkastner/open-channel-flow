% 2015-02-13 21:51:33.704599238 +0100
% Karl Kastner, Berlin
%% Dynamic Power Law Rating curve
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
classdef DynamicPowerRC < RatingCurve
	properties
		dh_dt0
		opt
		param0
		cext
		cflag
		np = 3;
	end % properties
	methods
		function obj = DynamicPowerRC(varargin)
			%obj.slope = slope;
			%obj.width = width;
			obj.cflag = false(obj.np,1);
			for idx=1:2:length(varargin)
				switch(varargin{idx})
				case {'afunc'}
					obj.afunc = varargin{idx+1};
				case {'pfunc'}
					obj.pfunc = varargin{idx+1};
				case {'c1'}
					obj.cflag(1) = true;
					obj.cext(1)  = varargin{idx+1};
				case {'c2'}
					obj.cflag(2) = true;
					obj.cext(2)  = varargin{idx+1};
				case {'c3'}
					obj.cflag(3) = true;
					obj.cext(3)  = varargin{idx+1};
				otherwise
					error('PowerRatingCurve','Unknown Parameter');
				end % switch
			end % for
			obj.rcfunc = @obj.rcfunc_;
			obj.opt = optimset( 'MaxIter', 10000, ...
				'MaxFunEvals', 10000, ...
				'Algorithm', 'levenberg-marquardt', ...
				'Jacobian', 'off', ... % on
				'TolFun', 1e-12, ...
				'TolX', 1e-12, ...
				'TolCon', 1e-12);
		end % PowerRatingCurve
		function [param obj] = fit_(obj,t0,q0,h0,dh_dt0)
			% lsqnonlin requires double input
			h0     = double(h0);
			q0     = double(q0);
			dh_dt0 = double(dh_dt0);

			% regression of exponential rating curve, start value is chosen from linear model
			% if no initial value is provided the exponential is fixed
			% to the theoretic value 1.5 and offset and scale are computed from
			% linear regression
			if (nargin() < 6)
				area = obj.afunc(h0);
				perimeter = obj.pfunc(h0);
				R = area./perimeter;

				% initial parameters from linearised fit
				A = [area.*R.^0.5];
				c = A \ q0;

				%A = [ones(size(h0)) h0];
				%c = A \ (q0.^(2/3));
				%param0 = [c(2).^(3/2) -c(1)/c(2) 1.5]';
				param0 = [c,c(1),0];
				param0 = param0(~obj.cflag);
				
			end
			%[param resn res flag] = lsqnonlin(@(c) W.*(obj.rcfunc(c,[],h0) - q0), param0, [],[],obj.opt)
			[param resn res flag] = lsqnonlin(@(c) obj.objective(c,t0,q0,h0,dh_dt0), param0, [],[],obj.opt)
		end % fit_

		function [Q gQ obj] = rcfunc_(obj,c_,t,h,dh_dt)
			%S  = obj.slope;
			%w  = obj.width;
			% (pseudo)-area
			A = double(obj.afunc(h));
			P = double(obj.pfunc(h));
			R = A./P;
			c = zeros(size(obj.cflag));
			% regressed, i.e. not externally provided parameters
			c(~obj.cflag,:) = c_;
			% externally provided parameters
			c(obj.cflag,:)  = repmat(cvec(obj.cext(obj.cflag)),size(c_,obj.np));

			% wave celerity
	
			% water velocity
			u = c(1)*R.^c(2);
			% stationary discharge
			Q = A.*u;
			% instationary discharge
			p   = c(3).*u.*dh_dt;
			fdx = (1 < p);
			p(fdx) = NaN;
			Q = Q.*sqrt(1 - p);
			

			% gradient of the discharge with respect to the parameters
			if (nargout() > 1)
				% TODO this is lazy, use analytic expression
				gQ = grad(@(c) obj.rcfunc_(c, t, h, dh_dt), cvec(c_));
			end % if nargout > 1
		end % rcfunc_
	end % methods
end % classdef DynamicPowerRatingCurve

