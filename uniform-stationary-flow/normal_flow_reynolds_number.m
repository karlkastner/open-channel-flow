% 2025-05-26 17:33:14.286205429 +0200
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
% function Re = normal_flow_reynolds(h,S,C)
function Re = normal_flow_reynolds_number(h,S,C)
	rho = 1000;
	g   = 9.81;
	mu  = 1e-3;
	u   = C*sqrt(h*S);
	Re  = rho*u*h/mu;
end

