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
% function q = flow_combined_resistance_discharge(h,S,c1,c2) %C,lcd)
function q = flow_combined_resistance_discharge(h,S,c1,c2) %C,lcd)
	if (~issym(h) && ~isa(h,'symfun'))
		g = Constant.gravity;
	else
		syms g positive;
	end
	if (issym(c2) || c2 ~= 0);
	den = sqrt(c1.*c1 + 4*g*c2.*abs(S).*h.*h.*h);
	% discharge at the interfaces
	% note that den > c1, so flow is always directed against the slope
	% except for the special case c2 = 0, which is captured below
	q  = sign(S).*(c1 - den)./(2*c2);
	else
		% cd1 q = g h S^2
		q=-g./c1*h.*h.*h.*S;
	end
end

