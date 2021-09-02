% Thu  3 Sep 15:10:39 +08 2020
% discharge on a suddenly tilted plate (suddenly unfrozen water on a tilted plate)
function Q = discharge_step_response(t,h,w,S,Cd)
	g = Physics.gravity;
	Q = -w.*sqrt(g.*h.^3.*S./Cd).*tanh(sqrt(S.*Cd.*g./h).*t);
%	if (isscalar(S) && 0 == S)
%%		Q(:) = 0;
%	else
%		Q(0 == S) = 0;
%	end
	if (isscalar(Cd) && 0 == Cd && S == 0)
		Q(:) = 0;
	else
		Q(0 == Cd & S == 0) = 0;
	end
end
