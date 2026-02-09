% 2024-09-27 15:22:58.952603374
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
function dh_dt = dh_dt(obj,t,h)
	if (obj.input_factor ~= 1)
	% TODO move input factor to RK
	h  = obj.input_factor * h;
	end
	if (0)
	dx = obj.dx;

	qi = obj.interface_values(t,h,0,0);

	% dh/dt = precipitation - infiltration - dqx/dx
	dh_dt = (qi(1:end-1) - qi(2:end))/dx;
	end

	% dh/dt = precipitation - infiltration - dqx/dx
	dh_dt = obj.dq_dx(t,h);
	dh_dt = cvec(dh_dt);

	if (~obj.opt.step_split_source)
		% this is for the entire time step, not sub-step for staged solvers
		dh_dt = dh_dt + obj.precipitation_rate(obj.aux.told,obj.aux.dt);

		ri = obj.infiltration_rate(h,dh_dt);
		dh_dt = dh_dt - ri;
	end

	if (obj.output_factor ~= 1)
	dh_dt = obj.output_factor*dh_dt/obj.input_factor;
	end
%	if (nargout>1)
%		qi = obj.output_factor*qi/obj.input_factor;
%	end
end

