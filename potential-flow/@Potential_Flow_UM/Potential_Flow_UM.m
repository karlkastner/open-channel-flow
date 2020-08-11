% Wed 20 May 11:40:20 +08 2020
%
%% numerical solver for flow on an unstructured mesh (triangulation)
%% by means of the Finite Element Method
%
% notation : small    s,n : flow
%	     captital S,N : mesh
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
classdef Potential_Flow_UM < Potential_Flow
	properties
		mode = 'lagrange';
	end % properties
%	methods (Access = protected)
%	end
	methods
		% constructor
		function obj = Potential_Flow_UM()
			obj.mesh = UnstructuredMesh();
			%obj.bcfun = @bcfun_dummy;
			%function [p,rhs] = bcfun_dummy(X,Y,id,val)
			%	p   = 0;
			%	rhs = 0;
			%end
		end % constructor
	end % methods
end % classdef Potential_Flow_UM

