% Wed  9 Aug 15:08:09 CEST 2017
%% process the transverse velocity profile
function [f] = process_transverse_profile(f,u,Q,A,x)
	neff = f.neff;

	f.t.scale = neff/(neff-(f.t.order+1));

	% weight
	f.t.weight = bsxfun(@times,u,rvec(Q).*rvec(A));
	w(:,1,:) = f.t.weight';

	% instantaneous profiles
	f.t.meas = bsxfun(@times,u,rvec(A)./rvec(Q));
	f.t.meas(:,0 == Q) = 0;

	% pre-filter
	f.t.meas_f = f.t.meas;
	if (f.nf.pre(1)>0)
		win     = kaiserwin(1:f.nf.pre(1));
		f.t.meas_f = wmedfilt(win,f.t.meas_f,1);
	end
	if (f.nf.pre(2)>0)
		win     = kaiserwin(1:f.nf.pre(2));
		f.t.meas_f = wmedfilt(win,f.t.meas_f',0)';
	end

	% weighted mean profile
	f.t.mean = wmean(f.t.weight,f.t.meas,2);

	% define inner region
	if (isempty(f.inner))
		f.inner = f.t.mean > 1;
	end

	% fit transverse profile
	f.t.poly = PolyOLS(f.t.order);
	f.t.poly.fit(x,f.t.meas_f',w);
	f.t.pred = f.t.poly.predict(x)';

	% prediction residual
	f.t.res = f.t.pred - f.t.meas;
	%f.t.res  = bsxfun(@minus,f.t.mean,f.t.A);
	
	% weighed rmse in profile centre
	f.t.s      = wrms(flat(f.t.weight(f.inner,:)),flat(f.t.res(f.inner,:)));

	% weighted rmse along profile
	% same as sqrt(scale*wvar(weight',f.t.A',1)');
        f.t.rmse = sqrt(f.t.scale)*wrms(f.t.weight,f.t.res,2);

	% pos-filter
	if (f.nf.post>0)
		win     = kaiserwin(1:f.nf.post);
		f.t.rmse = sqrt(wmedfilt(win,f.t.rmse.^2,1));
	end	

	% determine prediction error along range
	f.t = f.profile_prediction_error(f.t,true,Q);
end % process_transverse

