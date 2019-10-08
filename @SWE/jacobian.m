% Thu Apr 28 04:01:50 MSD 2011
% Karl Kästner
%
%% Jacobian of the SWE
%%
%% dq/dt + J dq/dx = sourceterm
%% note: d/dx(A*q) = J dq/dx
function [J, obj] = jacobian(obj,q)
        n = length(q)/2;
        h = q(1:n);
        m = q(n+1:end);
	Z = sparse(n,n);
	I = speye(n);
	if (~issym(q))
	g = obj.g;  
        J = [                             Z,            I;
             diag(sparse(-(m./h).^2 + g*h)), diag(sparse(2*m./h)) ];
	else
	syms g
        J = [                       Z,            I;
             diag((-(m./h).^2 + g*h)), diag((2*m./h)) ];
	end
	% note, used to be transposed
%	error ! see sw-matrix !

	% if QA
	% 
	%
end % function jacobian

