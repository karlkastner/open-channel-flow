% 2025-12-08 17:29:28.069899479
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
% TODO this is actually - dq/dx - dq/dy
function dh_dt = dq_dx(obj,t,h)
	[qxi,qyj] = obj.interface_values(t,h,false,false,false);

	dx  = obj.dx;

	% dh/dt = precipitation - infiltration - dqx/dx - dqy/dy
	dh_dt = ( (qxi(1:end-1,:) - qxi(2:end,:))/dx(1) ...
		+ (qyj(:,1:end-1) - qyj(:,2:end))/dx(2) ...
		);
end

