% Fri  9 Mar 17:56:18 CET 2018
%% side weir, along channel surface gradient
% x: distance from upstream end of side weir
function dzs_dx = dzs_dx(obj,x,z)
	g  = obj.g;
	alpha = obj.alpha;

	Qu = obj.Q.upstream;
	hu = obj.h.upstream;
	wu = obj.width.upstream;
	Qs = obj.Q.side;
	ws = obj.width.side;
	Au = hu*wu;

	% specific discharge of side channel
	switch (obj.shape)
	case {'rectangular'}
		dQ_dx = -Qs/ws;
	case {'triangular'}
		% TODO scale with h^3/2 than h^1
		dQ_dx = -Qs*2*x/ws^2;
	end
	% remaining discharge in main channel at x
	Q     =  Qu - (x/ws)*Qs;
%	Q = Q0;
	%dzs_dx = -alpha/g*(Q/A^2)*dQ_dx ...
	%	/ (1 - alpha/g*(Q^2*w0/A^3));
	% 12.8 and 12.16 in chow
	% A.157 in rosier 2007
	% TODO there seems to be a factor of two missing
	dzs_dx = -alpha*(Q*hu)*dQ_dx ...
		/ (g*hu.*Au^2 - alpha*Q^2);
end

