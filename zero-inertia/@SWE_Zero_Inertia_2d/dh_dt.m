% Thu  3 Oct 17:14:25 CEST 2024
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
function [dh_dt, qxi, qyj] = dh_dt(obj,t,h)
	% TODO move input factor to RK
	h = obj.input_factor * h;

	h = reshape(h,obj.nx);

	if (0)
	[qxi,qyj] = obj.interface_values(t,h,false,false,false);

	% dh/dt = precipitation - infiltration - dqx/dx - dqy/dy
	dh_dt = ( (qxi(1:end-1,:) - qxi(2:end,:))/dx(1) ...
		+ (qyj(:,1:end-1) - qyj(:,2:end))/dx(2) ...
		);
	end
	dh_dt = obj.dq_dx(t,h);

	if (~obj.opt.step_split_source)
	% this is for the entire time step, not sub-step for staged solvers
	dh_dt = dh_dt + obj.precipitation_rate(obj.aux.t, obj.aux.dt);

	ri    = obj.infiltration_rate(h,dh_dt);
	dh_dt = dh_dt - ri;
	end

	% flatten the vector
	dh_dt = dh_dt(:);
	dh_dt = obj.output_factor*dh_dt/obj.input_factor;
end % dh_dt

