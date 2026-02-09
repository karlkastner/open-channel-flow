% Wed 25 Mar 10:05:49 +08 2020
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
% c.f. engmann 1986
function cdq = darcy_weisbach2drag(f)
	if (issym(f))
		syms g positive
	else
		g = Constant.gravity;
	end
	% TODO this was 2g/f, why?
	C = sqrt(8*g./f);
	cdq = chezy2drag(C);
end
