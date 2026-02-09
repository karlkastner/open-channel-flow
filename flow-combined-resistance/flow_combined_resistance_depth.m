% 2025-05-13 12:02:00.554396228 +0200
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
% function h = flow_combined_resistance_depth(q,S,c1,c2)
function h = flow_combined_resistance_depth(q,S,c1,c2)
	g = Constant.gravity;
	%h = (C*lcd - (g/C*u/(0.5*sign(S))).^2 - C*C*lcd*lcd)/(4*abs(S)*g^2);
	if (0)
		h = (g*u.^2 - lcd*sign(S)*C^2*u)/(C^2*g*abs(S));
	else
		h = cbrt((c2.*q.*q + c1.*abs(q))./(abs(S)*g));
	end
end

