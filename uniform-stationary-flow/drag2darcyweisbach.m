% 2017-09-14 13:51:14.099671825 +0200
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
%% convert drag coefficient to chezy coefficient
%% g dz_s/dx + cd w u^2/h   = 0 (swe formalism)
%%      - S  +  1/C^2 U^2/H = 0 (chezy formalism)
function f = drag2darcyweisbach(cdq)
	if (issym(cd))
		syms g
	else
		g = Constant.gravity;
	end
	C = sqrt(g./cdq);
	f = chezy2darcy_weisbach(C);
end

