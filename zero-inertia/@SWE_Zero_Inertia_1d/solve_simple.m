% 2025-04-08 12:56:21.082462194 +0200
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
function [hh,exitflag,rmsr,rmsg,res,iter,out,iter_C] = solve_simple(obj,h0,T)
	hh = zeros(numel(h0),length(T));
	hh(:,1) = flat(h0);
	exitflag = zeros(length(T),1);
	% output counter
	obj.aux.odx = 1;
	obj.aux.timer = tic();
	for idx=2:length(T)
%[idx,length(T)]
		obj.aux.tdx=idx;
		dt = T(idx)-T(idx-1);
		%[rmsr,rmsg,res,out{idx},iter_C{idx}] = obj.step_simple(T(idx-1),hh(:,idx-1),dt);
		[hh(:,idx),stat] = obj.step_simple(T(idx-1),hh(:,idx-1),dt);

		exitflag(idx,1) = stat.flag;
		iter(idx,:) = stat.iter;
		iter_C{idx} = {stat.linear_solver.iter};
		rmsg = stat.rmse;
		rmsr = stat.rmse;
		res  = [];
		out{idx} = stat;

		if (exitflag(idx))
			hh = hh(:,1:idx);
			exitflag = exitflag(1:idx);
			disp(sprintf('no convergence at time %g',T(idx)))
			obj.finish(hh(:,end));			
			break;
		end
	end
end % solve_simple

