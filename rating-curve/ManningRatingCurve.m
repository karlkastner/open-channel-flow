% Mon Mar 16 09:59:33 CET 2015
% Karl Kastner, Berlin
%
% Stationary Rating Curve, Manning roughness formulation
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
classdef ManningRatingCurve < RatingCurve
	properties
		opt
		param0
		cext
		cflag
	end % properties
	methods
		function obj = ManningRatingCurve(varargin)
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
		end % constructor
		function [param obj] = fit_(obj,t0,q0,h0, dh_dt0)
			% The Manning formula has two parameters
			% energy slope and roughness, however, these are
			% are linearly dependent and combine to just one effective parameter
		
			area      = obj.afunc(h0);
			perimeter = obj.pfunc(h0);
			radius    = area./perimeter;

			% gradient (regression matrix)
			A         = [area.*radius.^(2/3)];
			param     = A \ q0;
		end % fit_

		function [Q G obj] = rcfunc_(obj, c, t, h, dh_dt)
			area       = obj.afunc(h);
			perimeter  = obj.pfunc(h);
			radius     = area./perimeter;

			% gradient (prediction matrix)
			G         = [area.*radius.^(2/3)];

			% discharge
			Q = G*c;

		end % rcfunc_
	end % methods
end % classdef ManningRatingCurve

