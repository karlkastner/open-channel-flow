% 2015-02-13 21:51:33.704599238 +0100
% Karl Kastner, Berlin
%
%% stationary rating curve, power law
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
classdef PowerRatingCurve < RatingCurve
	properties
		opt
		param0
		cext
		cflag
		np = 2;

	end % properties
	methods
		function obj = PowerRatingCurve(varargin)
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
		function [param obj] = fit_(obj,t0,q0,h0)
			% lsqnonlin requires double input
			h0 = double(h0);
			q0 = double(q0);

			% regression of exponential rating curve, start value is chosen from linear model
			% if no initial value is provided the exponential is fixed
			% to the theoretic value 1.5 and offset and scale are computed from
			% linear regression
			if (nargin() < 5)
				%A = [ones(size(h0)) h0];
				%c = A \ (q0.^(2/3));
				%param0 = [c(2).^(3/2) -c(1)/c(2) 1.5]';
				% initial parameters are determined by setting the exponent to 1
				area      = obj.afunc(h0);
				perimeter = obj.pfunc(h0);
				R         = area./perimeter;
				%c    = [area.*h0 area] \ q0;
				%param0 = [c(1);-c(2)/c(1);1];
				% fix offset to zero and exponent to 0.5
				% TODO linearised fit (log-space) for initial value
				A      = area.*R.^(0.5);
				c      = A \ q0;
				param0 = [c; 0.5];
				param0 = param0(~obj.cflag);
			end
			param0 = double(param0);
			[param resn res flag] = lsqnonlin(@(c) obj.objective(c,t0,q0,h0), param0, zeros(size(param0)),[],obj.opt);
		end % fit_

		function [Q gQ obj] = rcfunc_(obj, c_, t, h)
			% (pseudo)-area
			A = double(obj.afunc(h));
			P = double(obj.pfunc(h));
			R = A./P;
			c = zeros(size(obj.cflag));
			% regressed, i.e. not externally provided parameters
			c(~obj.cflag,:) = c_;
			% externally provided parameters
			c(obj.cflag,:)  = repmat(cvec(obj.cext(obj.cflag)),size(c_,obj.np));

			% discharge
			%f = (c(1,:).*A.*abs(R-c(2,:)).^c(3,:)).*sign(R - c(2,:));
			Q = c(1,:).*A.*R.^c(2,:);

			% gradient of the discharge with respect to the parameters
			if (nargout() > 1)
				% TODO this is lazy, use analytic expression
				gQ = grad(@(c) obj.rcfunc_(c, t, h), cvec(c_));
%				d   = cvec(R)-c(2);
%				sig = sign(d);
%				d   = abs(d);
%				g   = [  sig.*A.*d.^c(3), ...
 %                                       -c(1)*c(3)*sig.*A.*d.^(c(3)-1), ... %c(1)*c(3)*sign(d).*abs(d).^(c(3)-1), ...
  %                                       c(1)*sig.*A.*log(d).*d.^c(3) ];
%				g = g(:,~obj.cflag);
			end % nargout
		end % rcfunc_
	end % methods
end % classdef PowerRatingCurve

