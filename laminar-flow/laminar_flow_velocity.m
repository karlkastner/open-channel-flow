% 2022-04-06 18:11:09.054664099 +0200
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
%
function u = laminar_flow_velocity(h,S,c1)
	% g h^2 S + c1 u = 0
	% u = -g h^2 S/c1
	if (~issym(h))
		g  = 9.81;
		nu = Constant.viscosity_kinematic_water;
	end
	%if (nargin()<3)
	%end
	%u  = g/(3*nu)*S.*h.^2;
	u = -g/c1*h.^2.*S;
end
 


