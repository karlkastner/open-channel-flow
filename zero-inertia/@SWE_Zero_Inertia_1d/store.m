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
function store(obj,t,z)
	% call superclass storage function
	store@RAD_Model(obj,t,z);

	odx = obj.aux.odx;
	no = size(obj.out.h,1); %length(obj.out.esum);
	if (obj.opt.output.store_fluxes)
	if (odx > no)
		% reallocate
		no = 2*no;
		obj.out.h(no,1) = 0;
		obj.out.flow_x(no,1) = 0;
		obj.out.flow(no,1) = 0;
		obj.out.precipitation(no,1) = 0;
		obj.out.duration(no,1) = 0;
		obj.out.n_step_swe(no,1) = 0;
		obj.out.n_step_zie(no,1) = 0;
		obj.out.twet(no,1) = 0;
		obj.out.infiltration(no,1) = 0;
		obj.out.celerity_x(no,1) = 0;
		obj.out.diffusion_x(no,1) = 0;
		% TODO normalize diffusion of last time step by dt
	end
	end
end % store

