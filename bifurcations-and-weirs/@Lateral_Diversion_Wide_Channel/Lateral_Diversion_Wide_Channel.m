% Thu 30 Aug 18:12:46 CEST 2018
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
classdef Lateral_Diversion_Wide_Channel < Potential_Flow_Analytic
	properties
		% scaled velocity of diverted flow : v0 = alpha*uin
		alpha

		% scaled depth : beta ~ fs h/ws
		beta

		% secondary flow strength factor
		fs = 2/Constant.Karman^2;

		% scale of x-coordinate
		ws

		% scale of flow velocity
		uin

		shape = 'const';
		funfilename  = 'ldwc-functions.mat';
	end % properties
	
	methods
		function obj = Lateral_Diversion_Wide_Channel(varargin)
			obj = obj@Potential_Flow_Analytic();
			for idx=1:2:length(varargin)
				obj.(varargin{idx}) = varargin{idx+1};
			end
		end

		function Qs = Qs(obj)
			Qs = obj.uin*obj.alpha*obj.h*obj.ws;
		end

		function h = h(obj)
			 h = obj.beta*obj.ws/obj.fs;
		end

		function scales(uin,h,Qs,ws)
			
			
		end

	end % methods
end % Lateral_Diversion_Finite_Width

