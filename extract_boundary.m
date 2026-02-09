% 2024-11-02 09:31:13.453303681 +0100
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
% ij2id : spatial index to patch index
% boundary.ij = i and j coordinate of boundary points
% boundary points are as of now not yet chained to a line
function boundary = extract_boundary(ij2id)
	n = size(ij2id);
	boundary.ij  = [];
	by = [];
	k = 0;
	for idx=1:n(1)-1
		for jdx=1:n(2)-2
			if (ij2id(idx,jdx) ~= ij2id(idx+1,jdx))
				%xb(idx,jdx) = 1;
				k = k+1;
				boundary.ij(k,1:2) = [idx+1/2,jdx];
				%boundary.cid     = [ij2id(idx,jdx)
			end
			if (ij2id(idx,jdx) ~= ij2id(idx,jdx+1))
				k = k+1;
				boundary.ij(k,1:2) = [idx,jdx+1/2];
				%yb(idx,jdx) = 1;
			end
		end % for idx2
	end % for id1
end % extract_boundary

