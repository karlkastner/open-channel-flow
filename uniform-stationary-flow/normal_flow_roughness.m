% Sun 23 May 10:27:07 CEST 2021
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
% Karl Kastner, Berlin
%
%% roughness coefficient from uniform stationary flow
%function rgh = normal_flow_roughness(Q,H,W,S,type)
function rgh = normal_flow_roughness(Q,H,W,S,type)
	U = Q./(H.*W);
	rgh = U./sqrt(H.*S);	

	if (nargin()>4)
	switch (lower(type))
	case {'n','manning'}
		rgh = chezy2drag(rgh,H);
	case {'chezy','cz'}
		% nothing to do
	case {'drag','cd'}
		rgh = chezy2drag(rgh);
	case {'darcy-weisbach'}
		rgh = chezy2darcy_weisbach(rgh);
	otherwise
		error('here');
	end
	end
end

