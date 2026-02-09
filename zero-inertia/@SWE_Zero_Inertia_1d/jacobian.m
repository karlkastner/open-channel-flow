% 2024-09-27 10:16:32.799454273 +0200
% Fri  4 Oct 11:41:53 CEST 2024
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
% J = d/dh(dh/dt)
function [J] = jacobian(obj,t,h,dummy)
	h    = obj.input_factor*h;
	flag = 0;
	dx   = obj.dx;
	h    = rvec(h);
	h    = real(h);
	zb   = rvec(obj.zb);
	g    = obj.g;

	if (obj.opt.analytic_jacobian)
		obj.interface_values(t,h,false,true);

		J = -[ -obj.aux.dqxi_dhl(1:end-1); 
                        (obj.aux.dqxi_dhl(2:end) - obj.aux.dqxi_dhc(1:end-1));
                        obj.aux.dqxi_dhc(2:end)]/dx;
		J(:,obj.aux.ri_limited) = 0;
	
		if (~isempty(obj.infiltration_rate_linear))
			J(2,:) = J(2,:) + rvec(obj.infiltration_rate_linear);
		end

		J = reshape(J',obj.nx,3);
		if (obj.opt.returnmat)
			% TODO store precomputed indices
			J = diags2mat1d(t,J,obj.nx,obj.boundary_condition{1},obj.boundary_condition{2},dx);
		end % if returnmat
		J = obj.output_factor*J;
	else
		J = jacobian_regular_grid_1d(@(t,h) cvec(obj.dh_dt(t,cvec(h))),t,h,obj.nx,obj.opt.Jeps,obj.opt.twosided);
	end
	% this is not necessary in 1D
	obj.aux.max_eig_J = 0; %2*abs(min(J(5,:)));
end % function jacobian

