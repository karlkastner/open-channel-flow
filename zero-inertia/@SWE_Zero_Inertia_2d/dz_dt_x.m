% Thu  3 Oct 17:14:25 CEST 2024
% Mon  7 Oct 21:12:44 CEST 2024
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
% note: the function returns dh/dt, is is named dz_dt to be comptible
%       child of the superclass RAD_Model
function dh_dt = dz_dt_x(obj,t,h,p)
	% TODO move input factor to RK
	h = obj.input_factor * h;
	dx  = obj.dx;

	h = reshape(h,obj.nx);

	qxi = obj.interface_values_x(t,h,false);

	dh_dt = (qxi(1:end-1,:) - qxi(2:end,:))/dx(1);

	% this is for the entire time step, not sub-step for staged solvers
	dh_dt = dh_dt + obj.precipitation_rate(obj.aux.t, obj.aux.dt);

	ri    = obj.infiltration_rate(h,dh_dt);
	dh_dt = dh_dt - ri;

	% flatten the vector
	dh_dt = dh_dt(:);
	dh_dt = obj.output_factor*dh_dt/obj.input_factor;
end % dh_dt_i

