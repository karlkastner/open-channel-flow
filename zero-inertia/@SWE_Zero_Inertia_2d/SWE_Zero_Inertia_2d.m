% 2025-05-15 12:44:39.151127035 +0200
% Karl Kastner, Berlin
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
classdef SWE_Zero_Inertia_2d < SWE_Zero_Inertia_1d
	properties
	end	
	methods
		function obj = SWE_Zero_Inertia_2d(varargin)
			obj = obj@SWE_Zero_Inertia_1d();
			obj.opt.output.base_str = 'zi2d-';
			obj.boundary_condition = {[0,0,1,0], ...
				      [0,0,1,0], ...
				      [0,0,1,0], ...
				      [0,0,1,0] ...
					};
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end
			end
		end
	function xx = xx(obj)
		dx=obj.dx;
		x = (0:obj.n(1)-1)'*dx(1);
		xx = repmat(x,1,obj.n(2));
	end
	end
end

