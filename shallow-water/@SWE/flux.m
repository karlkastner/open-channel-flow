% Thu Apr 28 03:59:41 MSD 2011
% Sat Feb 19 19:59:39 MSK 2011
% Karl Kästner, KTH Stockholm
%
%% st venant's shallow water equation fluw
% SWE::flux
%
% TODO, this seems not correct when [A,Q] instead of [h,q] (unit width) are used
function [f, A, obj] = flux(obj, t, q)
	g = obj.g;
	n = length(q)/2;
        h = q(1:n);
        m = q(n+1:end);
        % TODO limit h, to avoid div by zero
	if (nargout() > 1)
		if (~issym(q))
			Z = sparse(n,n);
			I = speye(n,n);
			sp = @sparse;
		else
			Z = 0;
			I = 1;
			sp = @(x) x;
		end
		% TODO no magic numbers
		if (1)
		% variant D A q	
			% A = [ 0, I;
			%       
		      	A = [                Z,              I;
			     diag(sp(0.5*g*h)), diag(sp(m./h))];
		else
		% variant A D q
			% A = [               0,     I;
                        %      -((Q/h)^2 + g*h), 2*Q/h]
			A = [ Z, I;
			      diag(sp(-(m./h).^2 + g*h)), diag(sp(2*m./h))];
		end
% what is faster?
%                      sparse(M.N,M.N,-(m./h).^2 + g*h,n,n,n), sparse(M.N,M.N,2*m./h,n,n,n) ];
%                      sparse(M.N,M.N,-(m./h).^2 + g*h,n,n,n), sparse(M.N,M.N,2*m./h,n,n,n) ];
        	%A = [ M.Z, M.I; sparse(M.N,M.N,         0.5*g*h,n,n,n),   sparse(M.N,M.N,m./h,n,n,n) ];
	        f = A*q;
	else
		f = [ m; 
                      (m.*m)./h + 0.5*g*(h.*h)];
	end
end % SWE::flux

