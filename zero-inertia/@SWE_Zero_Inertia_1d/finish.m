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
%
function finish(obj,z)
	%if (~obj.aux.surface_flow)
	%	z(end+1:end+prod(obj.nx)) = 0;
	%end
	finish@RAD_Model(obj,z);
	if (obj.opt.output.store_fluxes)
		no = obj.aux.odx;
		obj.out.h = obj.out.h(1:no,:);
		obj.out.flow_x = obj.out.flow_x(1:no,:);
		obj.out.infiltration = obj.out.infiltration(1:no,:);
		obj.out.precipitation = obj.out.precipitation(1:no,:);
		obj.out.duration = obj.out.duration(1:no,:);
		obj.out.n_step_swe = obj.out.n_step_swe(1:no,:);
		obj.out.n_step_zie = obj.out. n_step_zie(1:no,:);
		obj.out.twet       = obj.out.twet(1:no,:);
		obj.out.celerity_x = obj.out.celerity_x(1:no,:);
		obj.out.diffusion_x = obj.out.diffusion_x(1:no,:);
	end
end

