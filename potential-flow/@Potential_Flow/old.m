	% inflow
	if (0)
		% linearly varying across cs
		for idx=1:n(1)
			AA(id,:)  = 0;
			AA(id,id) = 1;
		end
		X = obj.mesh.X;
		Y = obj.mesh.Y;
		% Y = X';
		% Y = Y';

		x = X(:,1);
		y = Y(:,1);
		bcval                      = afun(x,y);
		bb(1:n(1))                 = bcval;

		% outflow
		for idx=(n(2)-1)*n(1):n(2)*n(1)
			A(id,:) = 0;
			A(id,id) = 1;
		end
		x = X(:,end);
		y = Y(:,end);
		outval                     = bfun(x,y);
		bb((n(2)-1)*n(1)+(1:n(1))) = outval;

	end

