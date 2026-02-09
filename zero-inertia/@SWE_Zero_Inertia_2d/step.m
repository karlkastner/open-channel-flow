% Tue  8 Apr 11:07:46 CEST 2025
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
%function [h_dt,exitflag,rmsr,rmsg,res,iter,out] = step(obj,t,h0,dt)
function [h_dt,exitflag,rmsr,rmsg,res,iter,out,iter_C] = step(obj,t,h0,dt)
	siz = size(h0);
	if (siz(2)~=1)
		h0 = h0(:);
	end
	switch (obj.integration_scheme)
	case {'aid'}
		dt_0p5 = 0.5*dt;
		dx = obj.dx;
		if (~isfield(obj.aux,'I') || isempty(obj.aux.I))
			obj.aux.I = speye(prod(obj.n));
		end
		% implicit along x
		dh_dt_j = obj.dh_dt_j(t,h0);
		resfixed = (h0 + (dt_0p5)*dh_dt_j);
		% TODO interate 1 - dt J into Jacobian
		jfuni = @(h) (obj.aux.I - dt_0p5*obj.jacobian_i(t+dt_0p5,h));
		% TODO select quasi-tri-diagonal solver
		[h_dt_0p5,exitflag,rmsr,rmsg,res,iter(1),out] = gauss_newton(...
					@resfuni, ...
					h0, ...
					obj.input_factor*obj.gn_abstol, ...
					obj.gn_maxiter, ...
					obj.linear_solver_name, ...
					jfuni, ...
					obj.linear_solver_tol, ...
					obj.linear_solver_maxiter ...
					);
		% implicit along y
		dh_dt_i = obj.dh_dt_i(t+dt_0p5,h0);
		resfixed = (h0 + (dt_0p5)*dh_dt_i);
		%i = obj.interface_values_i(h_dt_0p5);
		%resfixed = (h0 + (dt_0p5/dx(1))*flat(qi(1:end-1,:) - qi(2:end,:)));
		jfunj = @(h) (obj.aux.I - dt_0p5*obj.jacobian_j(t+dt,h));
		% TODO select quasi-tri-diagonal solver
		[h_dt,exitflag,rmsr,rmsg,res,iter(2),out,iter_C] = gauss_newton(...
					@resfunj, ...
					h_dt_0p5, ...
					obj.input_factor*obj.gn_abstol, ...
					obj.gn_maxiter, ...
					obj.linear_solver_name, ...
					jfunj, ...
					obj.linear_solver_tol, ...
					obj.linear_solver_maxiter ...
					);
	otherwise
		[h_dt,exitflag,rmsr,rmsg,res,iter,out,iter_C] = step@SWE_Zero_Inertia_1d(obj,t,h0,dt);
	end
	if (siz(2)~=1)
		h_dt = reshape(h_dt,siz);
	end
function res = resfuni(h)
	 %qi  = obj.interface_values_i(h); -> note, need output factor
	 %res = h - (dt_0p5/dx(1))*flat(qi(1:end-1,:) - qi(2:end,:)) - 0*resfixed;
	 dh_dt_i = obj.dh_dt_i(t,h);
	 res = h - (dt_0p5)*dh_dt_i - resfixed;
end

function res = resfunj(h)
	%qj = obj.interface_values_j(h);
	%res = h - (dt_0p5/dx(2))*flat(qj(:,1:end-1) - qj(:,2:end)) - resfixed;
	dh_dt_j = obj.dh_dt_j(t+dt_0p5,h);
	res = h - (dt_0p5)*dh_dt_j - resfixed;
end

end % function step

