% 2024-09-27 15:22:58.952603374
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
% note: this is actually -dx/dx
function dqdx = dqdx(obj,t,h)
	dx = obj.dx;

	qi = obj.interface_values(t,h,false,false);

	dqdx = (qi(1:end-1) - qi(2:end))/dx;
end

