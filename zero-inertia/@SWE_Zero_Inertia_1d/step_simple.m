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
function [z_dt,stat] = step_simple(obj,t,z,dt)
	obj.aux.t  = t;
	obj.aux.told  = t;
	obj.aux.dt = dt;
	if (length(obj.nx)>1)
		obj.aux.zold = reshape(z,obj.nx);
	else
		obj.aux.zold = z;
	end
	switch (obj.opt.nlsolver.name)
	case {'gauss-newton'}
		[z_dt,stat] = obj.aux.fstep(t,z,dt);
	otherwise
		%z_dt = step_fixed_point(t,z,dt);
		[z_dt,stat]= fixed_point_iteration(@afun_fp,z,obj.opt.nlsolver,obj.opt.linear_solver);
	end
	obj.aux.stat = stat;
%	obj.aux.stat          = struct();
%	obj.aux.stat.maxe     = 0;
	obj.aux.step.n_neg    = 0;
	obj.aux.step.n_iter   = 0;
	obj.aux.step.n_attempt= 0;
	obj.aux.step.n_solver_failed=0;
	obj.aux.step.n_error_tolerance_exceeded=0;
	%obj.aux.odx = obj.aux.odx +1;
	obj.step_postprocess(z_dt);
	obj.store(t+dt,z_dt);

function [AA,rhs] = afun_fp(h)
	% TODO inhomogeneous boundary
	% TODO A requires different bc for zb and h (!)
	% there is Ae and A
	switch (obj.opt.time_integration_scheme)
	case {'euler-implicit'}
		A = dt*obj.fixed_point_matrix(t,h);
		%A(1,:) = 0;
		%A(end,:) = 0;
		Azb = Azbfun(A);
		AA       = (obj.aux.I - A);
		rhs      = h0 + Azb;
	case {'midpoint'}
		A = dt*obj.fixed_point_matrix(t,0.5*(h+h0));
		%A(1,:) = 0;
		%A(end,:) = 0;
		Azb = Azbfun(A);
		AA = (obj.aux.I - 0.5*A);
		rhs = h0 + Azb + 0.5*A*h0;
	case {'trapezoidal'}
		A0 = dt*obj.fixed_point_matrix(t,h0);
		A1 = dt*obj.fixed_point_matrix(t,h);
		Am = 0.5*(A0+A1);
		Azb = Azbfun(Am);
		%A0(1,:) = 0;
		%A0(end,:) = 0;
		%A1(1,:) = 0;
		%A1(end,:) = 0;
		%Azb = 0.5*(A0+A1)*cvec(obj.zb(2:end-1));
		%Azb(1) = Azb(2);
		%Azb(end) = Azb(end-1);
		rhs = h0 + Azb + 0.5*A0*h0;
		AA = (obj.aux.I - 0.5*A1);
	otherwise
		disp(obj.opt.integration_scheme);
		error('scheme not implemented');
	end
	
	function Azb = Azbfun(A)	
		if (1==length(obj.nx))
			Azb      = A*cvec(obj.zb(2:end-1));
			%A(1,:)   = 0;
			%A(end,:) = 0;
			%Azb      = A*cvec(obj.zb(2:end-1));
			Azb(1)   = Azb(2);
			Azb(end) = Azb(end-1);
		else
			% TODO quick fix
			obj.aux.I  = speye(prod(obj.nx));
			Azb        = A*flat(obj.zb(2:end-1,2:end-1));
			Azb        = reshape(Azb,obj.nx);
			Azb(:,1)   = Azb(:,2);
			Azb(:,end) = Azb(:,end-1);
			Azb(1,:)   = Azb(2,:);
			Azb(end,:) = Azb(end-1,:);
			Azb        = flat(Azb);
		end
	end % Azbfun
end % afun_fp

end % step_simple

