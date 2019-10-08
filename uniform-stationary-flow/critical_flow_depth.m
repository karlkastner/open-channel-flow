% 2017-09-04 18:17:27.414808491 +0200
%
%% critical flow depth in uniform stationary flow
function hc = critical_flow_depth(Q,W)
	if (issym(Q))
		syms g
	else
		g = Constant.gravity;
	end
	hc = (1/g*Q.^2/W.^2).^(1/3);
end

