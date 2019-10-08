% Thu 13 Jul 11:19:51 CEST 2017
% Karl Kastner, Berlin
%% compute discharge
function Q   = csdischarge(q,zs,zb,dw)
	 D   = bsxfun(@minus,rvec(zs),cvec(zb));
	 fdx = (D>0);
	 q(~fdx) = 0;
	 Q = dw*sum(q).';
end

