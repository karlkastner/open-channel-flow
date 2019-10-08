% 2017-09-05 00:41:54.131459023 +0800
%% normal flow depth in uniform stationary flow
function hn = normal_flow_depth(Q,W,C,S0)
	hn = (Q/(W*C*sqrt(S0)))^(2/3)
end

