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
% c1 : linear friction coefficient
% function u = flow_combined_resistance_velocity(h,S,c1,c2)
function u = flow_combined_resistance_velocity(h,S,c1,c2)
	if (~issym(h))
		g = Constant.gravity;
	else
		syms g;
	end
	den = sqrt(c1.*c1 + 4*g*c2.*abs(S).*h.*h.*h);
	% velocity at interfaces
	% note that den > c1, so flow is always directed against the slope
	% except for the special case c2 = 0, which is captured below
	u  = sign(S).*(c1 - den)./(2*c2.*h);
	% capture special case
	if (~issym(c2))
		fdx = (c2 == 0);
		u(fdx) = -g*(h(fdx).*h(fdx).*S(fdx))./c1(fdx);
	end
end

