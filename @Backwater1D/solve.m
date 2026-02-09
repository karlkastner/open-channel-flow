% 2016-04-08 11:15:55.270978358 +0200
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
%
%% solve the gradually varied flow equation (backwater equation)
%% C : chezy
%% W : width
%% Q : discharge
%% S : bed slope
%% y0 : surface elevation at outflow
%
% dy/dx = 1./(1-F^2)*(S0 - S_f  - 2*beta*Q/(g*A^2)*I - Q^2/(g*A^2)*dbeta_dx)
% TODO h0 should better be z0 and zb
% TODO lateral inflow I = 0
% function [x, h, zs, dh_dx, dzs_dx] = solve(obj,Q0,Qt,chezy,width,zb,z0_downstream,X)
%
function [x, h, zs, dh_dx, dzs_dx] = solve(obj,Q0,Qt,chezy,width,zb,z0_downstream,X)
	obj.X  = X;
	obj.Q0 = Q0;
%	if (isempty(Qt))
%		Qt = 0;
%	end
	obj.Qt_ = Qt;
	obj.chezy_ = chezy;
	obj.width_ = width;
	obj.zb_    = zb;

if (0)
	[x, zs] = obj.solver(@(x,h) obj.dzs_dx(x,h), X, h0, opt);
	zb = x*S0-h0;
	h  = zs-zb;
else
	% initial depth
	h0 = z0_downstream-min(obj.zb(X(1))); %(1));

	% solve backwater ODE
	% the equation is solved for the flow depth, not surface elevation,
	% as it is more numerically stable
	% depth remains small, elevation not
	% make initial step relatove
	sopt = obj.sopt;
	if (isfield(sopt,'InitialStep'))
		sopt.InitialStep = sopt.InitialStep/abs(X(2)-X(1));
	end

%	[x, h] = obj.solver(@(xi,h) obj.dh_dxi(xi,h), [0, 1], h0, sopt);
%	x=X(1)+(X(2)-X(1))*x;
	[x, h] = obj.solver(@(xi,h) obj.dh_dx(xi,h), X, h0, sopt);
	%[x h] = solver(@(x,h) backwater_dh_dx(x,h,Q,C,S0,W,widechannel), X, h0, opt);

	if (abs(X(end)-x(end))>abs(X(2)-X(1))*sqrt(eps))
		error('ode solver failed');
	end

	if (nargout() > 2)
		zs  = obj.zb(x) + h;
%		plot(x,[zs,obj.zb(x)])
%pause
	end
	
end % else of if 0

	if (nargout() > 3)
		dh_dx = obj.dh_dx(x,h,Q0,chezy,S0,W);
	end

	if (nargout() > 4)
		dzs_dx = dh_dx + obj.dzb_dx(x);
	end
end % backwater1d::solve

