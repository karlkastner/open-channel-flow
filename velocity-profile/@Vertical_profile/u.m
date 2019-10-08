% Do 15. Okt 18:11:07 CEST 2015
% Sat  7 Jan 14:41:29 CET 2017
% Karl Kastner, Berlin
%
%% predict velocity along the vertical based on profile
%
function [u, obj] = u(obj,S,h)
	n = size(S);

	param = obj.param;

	if (isrow(h) & iscolumn(S))
		Z = S*h;
	else
		Z = S.*h;
	end

	if (~issym(S) && ~issym(h) && ~issym(param))
		% allocate memory
		u = zeros(n);
	else
		syms u
	end

	for idx=1:n(2)
		% prediction matrix
		A = obj.regmtx(Z(:,idx),h(idx));
		% predict
		u(:,idx) = A*param(:,idx);
	end % for idx
end % u

