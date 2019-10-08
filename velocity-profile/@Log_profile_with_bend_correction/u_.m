% Fri 21 Sep 09:36:46 CEST 2018
%% streamwise velocity
function u = u_(obj,S,h,c)
	n = size(S);

	if (isrow(h) & iscolumn(S))
		Z = S*h;
	else
		Z = S.*h;
	end

	if (~issym(S) && ~issym(h) && ~issym(c))
		% allocate memory
		u = zeros(n);
	else
		syms u
	end

	for idx=1:n(2)
	u(:,idx) = (  c(1,idx)*log(Z(:,idx)) ...
		    + c(2,idx) ...
		    + c(3,idx)*Z.*log(Z(:,idx)) ...
		    + c(2,idx)*c(3,idx)/c(1)*Z(:,idx));
	end % for idx
end % u

