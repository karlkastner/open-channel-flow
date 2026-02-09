% 2025-05-13 12:00:39.611392059 +0200
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
% function q = flow_combined_resistance_slope(q,h,c1,c2)
function S = flow_combined_resistance_slope(q,h,cd1,cd2)
	if (~issym(h))
		g = Constant.gravity;
	else
		syms g positive;
	end
	S = -(cd2*abs(q).*q + cd1*q)./(g*h.*h.*h);
end

