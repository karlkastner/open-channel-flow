% 2025-07-23 18:34:37.740271444 +0200
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
function finish(obj,z)
	%if (~obj.aux.surface_flow)
	%	z(end+1:end+prod(obj.nx)) = 0;
	%end
	finish@SWE_Zero_Inertia_1d(obj,z);
	if (obj.opt.output.store_fluxes)
		no = obj.aux.odx;
		obj.out.flow_y = obj.out.flow_y(1:no,:);
		obj.out.celerity_y = obj.out.celerity_y(1:no,:);
		obj.out.diffusion_y = obj.out.diffusion_y(1:no,:);
	end
end

