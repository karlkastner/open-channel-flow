% 2015-02-13 19:06:11.715207218 +0100
% Karl Kastner, Berlin
%
% statinary rating curve, polynomial expression
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
classdef PolyRatingCurve < RatingCurve
	properties
		order
	end % properties
	methods
		function obj = PolyRatingCurve(order)
			obj.order = order;
			obj.rcfunc = @obj.rcfunc_;
		end % constructor
		function [q A obj] = rcfunc_(obj,param,t,h)
			% prediction
			A   = vander_1d(h,obj.order);
			q   = A*param;
		end % rcfunc_
		% parameter estimate
		function [param C obj] = fit_(obj,t0,q0,h0) %,C,param0)
			% covariance or weight matrix
			A  = vander_1d(h0,obj.order);
			%if (nargin() < 4)
				param = A \ q0;
			%else
			%	param = (A'*C*A) \ (A'*(C*q0));
			%end
		end % fit_
	end % methods
end % classdef PolyRatingCurve

