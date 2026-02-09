% 2025-10-28 16:39:40.308760503 +0100
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
function [dqxi_dhl,dqxi_dhc] = dqxi_dh(obj,t,h0)
	[qi0] = obj.interface_values(t,h0);
	epsh = sqrt(eps);
	dqxi_dhl = zeros(obj.nx(1)+1,1);
	dqxi_dhc = zeros(obj.nx(1)+1,1);
	for idx=1:3
			id = (idx:3:obj.nx(1))';
			h = h0;
			h(id) = h(id) + epsh;
			[qi] = obj.interface_values(t,h);
			dqi_dh = (qi - qi0)/epsh;
			dqxi_dhc(id)   = dqi_dh(id);
			dqxi_dhl(id+1) = dqi_dh(id+1);
	end
end

