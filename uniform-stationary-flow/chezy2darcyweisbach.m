% Wed 25 Mar 10:03:31 +08 2020
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
% c.f. Engmann 1986
function f = chezy2darcy_weisbach(C)
	if (issym(C))
		syms g
	else
		g = Constant.gravity;
	end
	% note, this was 2g/C^2, why?
	f = 8*g./C.^2;
end


