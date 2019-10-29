% Fr 31. Jul 11:04:11 CEST 2015
% Karl Kastner, Berlin

%function vert = vertical_profile_parameter(N, ln_z0, h_, z_, l_, w_, W, E, nf_, nmax, mode, order)
function trans = transversal_profile_parameter(N, fn, h_, w_, W, E, nmax, nf_, order)

	% only consider central part of the cross section
	inner  = abs(N) < nmax;
	n      = sum(inner);

	% measurement error (variance between campaigns)
	% this removes
	% turbulent fluctuation, adcp noise, gps noise
	% and displacement of measurement up or downstream
	res = fn*W'*E - fn;
	
%	[acorr acov] = acf_man(res(inner,:)*cvec(w_),[],[],1);
	[acorr acov] = autocorr_man4(res(inner,:),[],[],1);
	acov = acov*cvec(w_);

	% average velocity profile	
	trans.f0 = fn*W';
	trans.camp.res  = res;
	trans.camp.acov = acov;
	
	% transversal profile prediction
	trans.poly      = PolyOLS(order,false);
	[param res]     = trans.poly.fit(h_,fn',cvec(w_));
	res             = res';
	serr            = trans.poly.serr';
	rms_            = rms(res');

	% relative error
	relres          = res./fn;
	ns = size(fn,2);
	% correct for finite sample size
	relrms          = rms(relres,2)*sqrt(nf_/(nf_-order-1));

	% autocorrelation and autocovariance function
	%[acorr acov]     = acf_man(res(inner,:)*cvec(w_),[],[],1);
	[acorr acov]     = autocorr_man4(res(inner,:),[],[],1);
	acov             = acov*cvec(w_);

	% weighed averag over campaigns
%	acov = acov*cvec(w_);
	trans.total.res  = res;
	trans.total.acov = acov;

	trans.model.acov = trans.total.acov - trans.camp.acov;
	trans.model.res  = trans.total.res  - trans.camp.res;
	trans.model.serr = serr;
	trans.model.rms  = rms_;
	trans.model.relrms  = relrms;

	% prediction error
	% spatial correlation of the prediction error
	a = trans.model.acov; a = a/a(1);
	rho = ar1delay(a,nf_,1)

	% the covariances at until lag nf_-1 are spoiled by the spatial filter,
	% therefore s2 is estimated from the covariances at lag nf_
	acorr_ = acfar1(rho,n);
	s21  = trans.model.acov(nf_+1)/acorr_(nf_+1);
%*rho^(-nf_);
%	s21  = trans.model.acov(nf_+1)*rho^(-nf_);
	% this is now the variance of the first element, but not the process variance
	n = size(res,1);
	[acorr0 acov0] = acfar1(rho,n,0);
	s2 = s21/acov0;
	s2 = s21;

	trans.model.rho = rho;
	trans.model.s2  = s2;

	fprintf('Transversal profile: s2: %f rho %f\n',s2,rho);
end

