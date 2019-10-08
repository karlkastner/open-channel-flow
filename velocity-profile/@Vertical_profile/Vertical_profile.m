% Sat  7 Jan 14:51:27 CET 2017
% Karl Kastner, Berlin
%
%% vertical profile of the streamwise velocity, superclass
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
classdef Vertical_profile < handle
	properties (Abstract)
	end % properties (Abstract)
	properties
		param
		S2
		serr
		% maximum relative distance above bottom for fit
		smax = 1;
		% for replacement with symbolic value
		karman = Constant.Karman;
	end % properties
	methods (Static)
		A = regmtx(z,h);
	end
	methods
		% constructor
		function obj = Vertical_profile(varargin)
			for idx=1:2:length(varargin)-1
				obj.(varargin{idx}) = varargin{idx+1};
			end
		%	if (nargin() > 0)
		%		obj.param = param;
		%	end
		%	if (nargin() > 1)
		%		obj.S2 = S2;
		%	end
		end % Vertical_profile

		% pseudo parameter
		function s2param = s2param(obj,fdx)
			if (isempty(obj.S2))
				if (nargin()<2)
					n = size(obj.param,2);
				else
				if (islogical(fdx))
					n = sum(fdx);
				else
					n = length(fdx);
				end
				end
				s2param = zeros(size(obj.param,1),n);
			else
			if (nargin() < 2)
				s2param = diag3(obj.S2);
			else
				s2param = diag3(obj.S2(:,:,fdx));
			end
			end
		end
	end % methods
end % classdef Vertical_profile


