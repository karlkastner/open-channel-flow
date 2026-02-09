% Mon 19 Feb 16:07:09 CET 2018
% Mon 21 Oct 17:42:14 +08 2019
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
%% convert chezy to manning
% function n = chezy2manning(C,R)
function n = chezy2manning(C,R)
	n = 1./C.*R.^(1/6);
	if (~issym(R))
		n(R<0) = NaN;
	end
end


