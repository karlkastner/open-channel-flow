% Sun 20 May 12:26:47 CEST 2018
%% Dynamic solution of the shallow water equation (depth average, 2D)
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
classdef SWE_2d < handle
	properties
		mesh
		bcfun
		h
		qx
		qy
		zb

		g = Constant.gravity();		
		C = 60;	
	end % properties
	methods
		function obj = SWE_2d()
			obj.mesh = SMesh();
		end % constructor
		function umag = umag(obj)
			umag = hypot(obj.ux,obj.qx);
		end
		function ux = ux(obj)
			ux = obj.qx./obj.h;
		end
		function uy = uy(obj)
			uy = obj.qy./obj.h;
		end
	end % methods
end % SWE_2d

