% Thu 15 Mar 15:46:43 CET 2018
%% cut a rectangle from the domain
%% TODO, this requires also an adaptation of the derivative matrices
%%       -> step over to semi-unstructured mesh
function obj = cut_rectangle(obj,id1,id2)
	n = obj.mesh.n;
	m = prod(n)*[1 1];
%fdx    = r > r_min & theta < theta0(1);
	%id1    = find(fdx(:,1));
	%id2    = find(fdx(end,:));
	%id     = find(fdx);

%	figure(300);
%x = obj.x;
%clf

	A = obj.A;
	b = obj.b;
	id = flat( (id1(1):id1(2))'*ones(1,id2(2)-id2(1)+1) ...
	           + n(1)*ones(id1(2)-id1(1)+1,1)*(id2(1)-1:id2(2)-1));
	
	A(id,:) = 0;
	A(sub2ind(m,id,id)) = 1;
	b(id) = 0;

	% non-flow accross cut-region boundaries
	% TODO corners
	if (id1(1) > 1)
		% one row below
		id                    = (id1(1)-1) + (id2(1)-1:id2(2)-1)*n(1);
		A(id,:)               = 0;
		A(sub2ind(m,id,id))   = -1;
		A(sub2ind(m,id,id-1)) =  1;
		b(id) = 0;
	end
	if (id1(2) < n(1))
		% TODO
	end
	if (id2(1) > 1)
		% right neighbour
		id                       = (id1(1):id1(2)) + (id2(1)-2)*n(1);
		A(id,:)                  = 0;
		A(sub2ind(m,id,id))      = -1;
		A(sub2ind(m,id,id-n(1))) =  1;
		b(id) = 0;
	end
	if (id2(2) < n(2))
		% left neighbour
		id                       = (id1(1):id1(2)) + id2(2)*n(1);
		A(id,:)                  = 0;
		A(sub2ind(m,id,id))      = -1;
		A(sub2ind(m,id,id+n(1))) =  1;
		b(id) = 0;
	end
	obj.A = A;
	obj.b = b;
%	b = reshape(b);
end
