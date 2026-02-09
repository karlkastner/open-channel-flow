% 2025-07-30 10:03:24.958030227 +0200
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
% function k = chezy2klin(cz,h,S)
function k = chezy2klin(cz,h,S)
	% g*h*S + cd*u^2  -> cd = g h S/u^2 -> u = sqrt(g/cd*h*S)                         
	% g*h*S + k*u/h     -> k  = g h^2 S/u = cd*h*u 
	% 			u = g h^2 S/k 
	% ->  k = cd*h*sqrt(g/cd*h*S) = h*sqrt(g*cd*h*S)
	cdrag = chezy2drag(cz);
	g = 9.81;
	k = h*sqrt(cdrag*g*h*S);
end

