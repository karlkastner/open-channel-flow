% 2025-08-06 19:13:15.897720793 +0200
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
function init_solve(obj)
	init_solve@RAD_Model(obj);

	if (obj.opt.step_split_source)
		obj.aux.fstep2 = obj.aux.fstep;
		obj.aux.fstep  = @obj.step_split_source;
	end
	obj.aux.max_eigJ = sqrt(eps);

	if (obj.opt.output.store_fluxes)
		if (isfield(obj.opt.output,'no'))
			no = obj.opt.output.no;
		else
			no           = length(obj.out.to);
		end
		output_class_str = func2str(obj.opt.output.class);
		% make single to avoid overflow
		nxi = obj.nx;
		nxi(1) = nxi(1)+1;
		obj.out.h             = zeros(no,prod(obj.nx),output_class_str);
		obj.out.flow_x        = zeros(no,prod(nxi),output_class_str);
		obj.out.infiltration  = zeros(no,prod(obj.nx),output_class_str);
		obj.out.precipitation = zeros(no,1,output_class_str);
		obj.out.duration      = zeros(no,1,'single');
		obj.out.n_step_swe      = zeros(no,1,'single');
		obj.out.n_step_zie      = zeros(no,1,'single');
		obj.out.twet          = zeros(no,prod(obj.nx),'single');
		obj.out.diffusion_x   = zeros(no,prod(nxi),output_class_str);
		obj.out.celerity_x    = zeros(no,prod(nxi),output_class_str);
		obj.out.n_step_zie    = zeros(no,1,'single');
		obj.out.n_step_swe    = zeros(no,1,'single');
	end
end % init_solve

