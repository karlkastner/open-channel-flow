% Tue 15 Aug 09:19:07 CEST 2017
function f = process_joint(f,Q)
	f.tv.rho        = wcorr(flat(f.v.weight(f.inner,:)),flat(f.v.res(f.inner,:)), ...
				  flat(f.t.weight(f.inner,:)),flat(f.t.res(f.inner,:)));
	f.tv.range_rmse = sqrt(f.t.range_rmse.^2 + f.v.range_rmse.^2)';
	f.tv.range_R2   = 1-mean((cvec(f.tv.range_rmse)*rvec(Q)).^2,2)/var(Q);
end

