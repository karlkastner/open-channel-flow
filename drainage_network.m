% Sat  2 Nov 09:15:32 CET 2024
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
%% function [ij2c,sink] = drainage_network(zs)
%%
%% input:
%% 	zs : ni x nj, water surface elevation of in case of DEM ground elevation
%% output:
%%	ij2c : ni x nj map with index of the catchment the pixel belongs too
%%	sink : nc x 1, i and j coordinate of catchment by id
%
%
function [ij2c,sink,next] = drainage_network(zs)
	% TODO no magic numbers
	% e = 1e-3;
	
	n    = size(zs);
	nn   = prod(n);
	% neighbour to which water flows
	next = zeros(n(1),n(2),2);
	% final point to which water flows
	sink = zeros(n(1),n(2),2);
	
	% di
	d1 = [-1, 0, +1, -1, 1, -1, 0, +1];
	% dj
	d2 = [-1, -1, -1, 0, 0, 1, 1, +1];
	% slope magnitude, one on straight lines and 0.7 along diagonals
	s = [1/sqrt(2), 1, 1/sqrt(2), 1, 1, 1/sqrt(2), 1, 1/sqrt(2)];
	
	% for row
	for id1=1:n(1)
		% dislplay progres
		disp(id1/n(1));
		% for each column
		for id2=1:n(2)
			% drain water from point recursively
			drain(id1,id2);
		end % for id2
	end % for idx1
	% catchment indices
	sink = [flat(sink(:,:,1)),flat(sink(:,:,2))];
	[sink,~,ij2c] = unique(sink,'rows');
	ij2c = reshape(ij2c,n);
% end of main body

% recursive formulation
function drain(id1,id2)
	% process pixel, if not jet processed
	if (0 == next(id1,id2,1))
		% max negative gradient
		gmin = 0;
		% start with assuming curreng pixel is a sink
		next(id1,id2,:) = [id1,id2];
		sink(id1,id2,:) = [id1,id2];
		% for each neighbour, test if water flows to neighbour
		for nd1=1:8
			id1_ = mod(id1+d1(nd1)-1,n(1))+1;
			id2_ = mod(id2+d2(nd1)-1,n(2))+1;
			% gradient to nieghbour nd1
			g = s(nd1)*(zs(id1_,id2_)-zs(id1,id2));
			if (g < gmin)
				% this neighbour is lower than all others before
				gmin=g;
				next(id1,id2,:) = [id1_,id2_];
			else % if g  < gmin
			% if equal, choose smallest index
			% TODO this is not quite right
			if (g == gmin && (id1_ + (id2_-1)*n(1) < id1 + (id2-1)*n(1)))
				%gmin = g;
				next(id1,id2,:) = [id1_,id2_];
			end % if equal
			end % if g < gmin
		end % for nd1 
		% drain recursively towards sink
		drain(next(id1,id2,1),next(id1,id2,2));
		% sink of this neighbour
		sink(id1,id2,:) = sink(next(id1,id2,1),next(id1,id2,2),:);
	end % if not yet processed
end % drain

%function merge()
%	s_ = 
%	s_ = unique(s_);
%	for idx=1:ns
%		for jdx=1:ns
%			% walk from j to j, if a point inbetween exceeds delta, do not merge
%			% if reached j, merge
%			d1 = i1-i1_;
%			d2 = i2-i2_;
%			% TODO inverted
%			for i_ = 0:d1
%				dzs_ = zs() - zs(i1,i2);
%				if (dzs > e)
%					break;
%				end
%			end
%			if (dzs < e)
%				% merge both sinks
%				s_(s_ == i12_) = i12;
%			end
%		end
%	end
%end

end % drainage_network

