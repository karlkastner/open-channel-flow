% Thu Apr 28 16:42:03 MSD 2011
% Thu Apr 28 04:28:25 MSD 2011
% Karl Kästner
%
%% eigenvalues und vectors of the swe
% SWE::matrix
function [Lambda, R, Rinv, obj] = fluxmateig(obj,y_bar,w_bar)
	g = obj.g;
	n = length(y_bar)/2;

	if (~issym(y_bar))
		sp = @sparse;
	else
%		syms g
		sp = @(x) x;
	end

	if (~obj.QAflag)
		h = y_bar(1:n);
		u = y_bar(n+1:end)./h;
	else
		A = y_bar(1:n);
		Q = y_bar(n+1:2*n);
		h = A./w_bar;
		u = Q./A;
	end
	c = sqrt(g*h);

	% 15.36 in randall
	l_1 = (u - c);
	l_2 = (u + c);

	I = speye(n);
	Z = sparse(n,n);

	% scales, so that eigenvectors are normalized
	flag = 0;
	if (flag)
		s      = diag(sparse(1./sqrt([1 + l_1.^2; 1+l_2.^2])));
		si     = diag(sparse(1./sqrt([1 + l_2.^2; 1+l_1.^2])));
		scale  = 1;
	else
		s=1;
		si=1;
		scale = 1./(2*c);
		scale = [scale;scale];
		scale = diag(sp(scale));
	end

	l_1 = diag(sp(l_1));
	l_2 = diag(sp(l_2));

	Lambda = [ l_1, Z;
                   Z, l_2];
	
	if (nargout() > 1)
		% 15.37 in randall (normalization missing in randall)
		R      = [     I,   I;
       	                     l_1, l_2]*s;
%		R      = R*s;
	end
	if (nargout() > 2)
		% 15.38 in randall
%		if (0)
		Rinv  = si*[l_2, -I;
       	                   -l_1,  I]*scale;
%		else
%			Rinv  = si*[l_2, -I;
 %      	                	   -l_1,  I];
%		end
	end
end % shallow_water_matrix

