% Mon 25 Aug 14:25:58 CEST 2025
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
function du_dh = flow_combined_resistance_du_dh(h,S,c1,c2)
	if (~issym(h))
		g = Constant.gravity;
	else
		syms g;
	end
	% velocity at interfaces
	den = sqrt(c1.^2 + 4*g*c2.*abs(S).*h.^3);
	u   = sign(S).*(c1 - den)./(2*c2.*h);
	
	du_dh = -u./h - 3*g*abs(S).*h./den;
end % du_dh

