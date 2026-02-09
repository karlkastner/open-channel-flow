% 2025-09-09 12:48:27.831974714 +0200
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
function step_postprocess(obj,hnew)
	told = obj.aux.told;
	dt = obj.aux.dt;
	% call parent class function
	step_postprocess@SWE_Zero_Inertia_1d(obj,hnew);
	if (obj.opt.output.store_fluxes)
		% TODO account for stages
		odx = obj.aux.odx;
		obj.out.flow_y(odx,:) = obj.out.flow_y(odx,:) + dt.*rvec(obj.aux.qyj(:));
		obj.out.celerity_y(odx,:) = obj.out.celerity_y(odx,:) + dt.*rvec(obj.aux.aj(:));
		obj.out.diffusion_y(odx,:) = obj.out.diffusion_y(odx,:) + dt.*rvec(obj.aux.dj(:));
	end
end

