%% Fri Feb 13 10:02:52 CET 2015
% Karl Kastner, Berlin
%% rating curve superclass
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
classdef RatingCurve < handle
	properties
		% nc,1 time at calibration
		t0
		% nc,1 stage at calibration
		h0
		% nc,1 discharge at calibraction
		q0
		% calibration parameters and goodness of fit with linearisation
		lin = struct('param',struct(),'res0',[],'serr0',[],'R2',[]);
		% goodness of fit for jackknife
		jk  = struct('param',struct(),'res0',[],'serr0',[],'R2',[]);
		% TODO replace rcfunc function pointer with a variable argument list
		% and rename rcfunc_ into rcfunc in derived classes
		rcfunc
		afunc = @(h) 1; % implicitely considered (width factors into c1 and depth increases c3 by 1)
		pfunc = @(h) 1;
	end % properties
	methods
		function [n obj] = nparam(obj)
			n = length(obj.lin.param);
		end
		function [n obj] = nsample(obj)
			n = length(obj.q0);
		end
		function obj = RatingCurve()
			% fit parameters with Jacknife
			obj.jk.param = Jackknife(@obj.fit_);
		end % constructor
		function obj = fit(obj,t0,q0,h0,varargin)
			obj.t0 = t0;
			obj.q0 = q0;
			obj.h0 = h0;
			obj.lin.param.val = obj.fit_(t0,q0,h0,varargin{:});
			% goodness of fit with linearisation
			obj.lingoodness(varargin{:});
		end
		function obj = jkfit(obj,t0,q0,h0,varargin)
			obj.t0 = t0;
			obj.q0 = q0;
			obj.h0 = h0;
			obj.jk.param.estimate(t0,q0,h0,varargin{:});

			% goodness of fit with jaccknife
			obj.jkgoodness(varargin{:}); %farg,jarg);
		end
		function [res g obj] = objective(obj, c, t0, q0, h0, varargin)
			%W = ones(length(h0),1);
			if (1 == nargout)
				[q]   = obj.rcfunc(c,t0,h0,varargin{:});
			else
				[q g] = obj.rcfunc(c,t0,h0,varargin{:});
			end
			%res = W.*(q - q0);
			res = q - q0;
		end
		function [q A bias serr sp obj] = predict(obj,t,h,varargin)
			if (1 == nargout())
				q = obj.rcfunc(obj.lin.param.val,t,h,varargin{:});
			else
				[q A] = obj.rcfunc(obj.lin.param.val,t,h,varargin{:});
			end % if
			bias = [];
			if (nargout() > 2)
				% 1-sigma confidence interval
				c = 1;
				t = tinv(normcdf(1),obj.nsample-obj.nparam);
				s2   = sum((A*obj.lin.param.C).*A,2);
				% error of mean response
				serr = t*sqrt(abs(s2));
				% predicition error
				sp   = t*sqrt(obj.lin.serr0.^2+abs(s2));
			end
		end % predict
		function [val A bias serr C obj] = jkpredict(obj,t,h,varargin)
			if (nargout() > 4)
				[val bias s2 C] = obj.jk.param.apply(obj.rcfunc,t,h,varargin{:});
			else
				[val bias s2] = obj.jk.param.apply(obj.rcfunc,t,h,varargin{:});
			end
			serr = sqrt(s2);
			A = [];
		end
		function obj = jkgoodness(obj,varargin)
		        [val void bias serr Cp] = obj.jkpredict(obj.t0,obj.h0,varargin{:});
%			res0      = val - obj.q0;
%			serr0     = sqrt(res0'*res0/(length(obj.q0)-length(obj.jk.param.val0)));
			% as this is a non-linear function the res for bias and serr differ
			res0  = [];
			serr0 = sqrt(mean(serr.^2));

			% coefficient of determination
			R2           = 1 - serr0*serr0/var(obj.q0);
			% write back
			obj.jk.res0  = res0;
			obj.jk.serr0 = serr0;
			obj.jk.R2    = R2;
			n            = length(obj.q0);
			k            = length(obj.jk.param.val0);
			obj.jk.aic   = akaike_information_criterion(serr0,n,k);
			obj.jk.bic   = bayesian_information_criterion(serr0,n,k);
		%	obj.jk.C     = C;
		end
		% goodness of fit
		function obj = lingoodness(obj,varargin)
			[val A]   = obj.predict(obj.t0,obj.h0,varargin{:});
			res0      = val - obj.q0;
			serr2     = sum(res0.*res0)/(obj.nsample-obj.nparam);
			serr0     = sqrt(serr2);
			% covariance matrix
			C = inv(A'*A);
			% error covariance matrix
			obj.lin.param.C   = serr2*C;
			% standard error of the parameters (not yet Studentised)
			obj.lin.param.serr = sqrt(diag(obj.lin.param.C));
			% coefficient of determination
			obj.lin.R2        = 1 - serr2/var(obj.q0);
			% hat matrix
 			obj.lin.param.hat = A*C*A';
			% write back
			obj.lin.res0  = res0;
			obj.lin.serr0 = serr0;
			obj.lin.aic   = akaike_information_criterion(serr0,obj.nsample,obj.nparam);
			obj.lin.bic   = bayesian_information_criterion(serr0,obj.nsample,obj.nparam);
		end % goodness
		function [ret D A] = testlingrad(obj,varargin)
			[val A]   = obj.predict(obj.t0,obj.h0,varargin{:});
			numA = grad(@(c) cvec(obj.rcfunc_(c,obj.t0,obj.h0)),obj.lin.param.val);
			%numA = grad(@(h) cvec(obj.predict(obj.t0,h,varargin{:})),obj.h0);
			D = (A-numA)./(A+numA)
			A
			if (max(abs(D(:))) < 10*sqrt(eps))
				disp('test passed')
				% POSIX_SUCCESS
				ret = 0;
			else
				disp('test failed')
				% posix error
				ret = 127;
			end
		end
	end % methods
end % classdef RatingCurve

