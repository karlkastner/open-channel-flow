% Tue 21 Apr 08:55:19 +08 2020
% discharge in a compound cross-section
% under uniform-stationary flow
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
classdef Compound_Cross_Section
	methods (Static)
		[A,h]   = area(zs,zb,width);
		[Q, qn] = discharge(zs,zb,W,C,S,ismanning);
		ret     = roughness(Q,zs,zb,W,S,ismanning);
	end
end

