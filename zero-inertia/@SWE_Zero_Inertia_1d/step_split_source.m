% 2025-10-28 19:04:16.635617665 +0100
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
% note that we ignore here for the time being the error introduced by varying infiltration and precipitation
function [h,stat] = step_split_source(obj,t,h0,dt)
	% ingerate source terms first half step
	ri = obj.infiltration_rate(h0,[]);
	rp = obj.precipitation_rate(t,0.5*dt);
	h  = max(0, h0 + 0.5*dt*(rp - ri));
	% step the flow entire step
	[h,stat] = obj.aux.fstep2(t,h,dt);
	% integrate source terms second half step
	ri = obj.infiltration_rate(h,[]);
	rp = obj.precipitation_rate(t+0.5*dt,0.5*dt);
	h  = max(0, h + 0.5*dt*(rp - ri));


	if (nargout()>1)
	% estimate the error
	tol = obj.time_integration_tolerance(h);
		
	% as constant rates are linear, they do not contribute to the error
	% this avoids computing flux dependent constant rates
	%obj.aux.ignore_constant_rates = true;
	dzdt_0 = obj.dz_dt(t,h0);
	dzdt_1 = obj.dz_dt(t+dt,h);
	%obj.aux.ignore_constant_rates = false;
	delta_dzdt     = (dzdt_1 - dzdt_0);
	[delta_dzdtmax,idmax] = max(abs(delta_dzdt));
	dzmax = dt*delta_dzdtmax;
	% TODO this constant is for trap, not for split
	C      = 0.12;
	emax   = C*dzmax;
	% safety margin is factored in later
	dt_opt = dt*sqrt(tol/emax);

	stat.rmse   = NaN;
	stat.maxe  = emax;
	stat.idmax = idmax;
	stat.dt_opt = dt_opt;
	end
end

