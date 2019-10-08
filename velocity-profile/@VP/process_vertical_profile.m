% Fr 31. Jul 11:04:11 CEST 2015
% Fri 11 Aug 10:28:31 CEST 2017
% Karl Kastner, Berlin
%
%% predict vertical profile error distribution parameter for HADCP error estimate
%

% ln_z0 : [nn  x nt] roughness length of individual transect
% zb    : [nn  x nc] bed elevation during each campaign
% zi    : [1]       instrument level (constant)
% zs    : [nc x  1] water level during each campaign
% inner : analysis is limited to the inner region of the cross section, where the data is well defined
% order : degree of polynomial to predict the profile (0==constant, 1==linear)
%
% note: an earlier version of this script could treat cross sections individual,
% and then subtracted the inter-calibration residual or variance from the total residual/variance,
% however, this yields similar results, so approach dropped for simplicity
%
function f = vertical_profile(f,us,ln_z0,pp,zb,zi,zs,Q,A,win,winrmse)
	neff = f.neff;

	f.v.scale = neff/(neff-(f.t.order+1));

	n = size(us);
	l = Log_profile_with_wake();
	l.shear_velocity(flat(us));
	l.ln_z0(flat(ln_z0));
	l.perturbation_parameter(flat(pp));

	% flow depth (h)
	h = bsxfun(@minus,rvec(zs),zb);

	% instrument depth
	hi = zi - zb;
	Si = bsxfun(@times,hi,1./h);
	
	% velocity at instrument depth (u)
	% TODO rename to u
	u  = l.u(flat(Si)',flat(h)')';

	% depth average velocity (bar u)
	ubar = l.ubar(flat(h));

	f.v.df_dln_z0 = real(reshape(l.df_dln_z0(flat(h),flat(hi)),n));

	f.v.u    = real(reshape(u,n));
	f.v.ubar = real(reshape(ubar,n));

	% average of the transverse velocity profile scale (u/bar u)
	f.v.meas    = f.v.u./f.v.ubar;

	f.v.weight = bsxfun(@times,f.v.u,rvec(Q).*rvec(A));
	w(:,1,:)   = f.v.weight';

	% average profile
	f.v.mean = wmean(f.v.weight,f.v.meas,2);

	x = zs;
	switch (lower(f.v.mode)) % prediction from roughness lenght and wake parameter
	case {'profile'}
		ln_z0_f = ln_z0;
		pp_f    = pp;
		% filter the roughness length spatially
		fprintf('Before filtering: std %f skew %f kurt %f\n', ...
			qstd(flat(ln_z0)), qskew(flat(ln_z0)), qkurtosis(flat(ln_z0)));

		if (f.nf.pre(1)>0)
			win     = kaiserwin(1:f.nf.pre(1));
			ln_z0_f = wmedfilt(win,ln_z0_f,1);
			pp_f    = wmedfilt(win,pp_f,1);
		end
		if (f.nf.pre(2)>0)
			win     = kaiserwin(1:f.nf.pre(2));
			ln_z0_f = wmedfilt(win,ln_z0_f',0)';
			pp_f    = wmedfilt(win,pp_f',0)';
		end
		% Note : the standard variation is not so much reduced, because
		% of the systematic differences between the campaigns
		fprintf('After filtering: std %f skew %f kurt %f\n', ...
		qstd(flat(ln_z0)), qskew(flat(ln_z0)), qkurtosis(flat(ln_z0)));

		%ln_z0 = wmean(f.v.weight,ln_z0,2)*ones(1,n(2));
		%pp    = wmean(f.v.weight,pp,2)*ones(1,n(2));
		%nf = [1 3];
		%u = medfilt1(medfilt1(u,nf(1))',nf(2))';
		%ubar = medfilt1(medfilt1(ubar,nf(1))',nf(2))';

		% shear velocity does not need to be predicted
		% predict log of roughness length for each campaign
		f.v.poly.ln_z0 = PolyOLS(f.v.order);
		f.v.poly.ln_z0.fit(x,ln_z0_f',w); 
		f.v.poly.pp = PolyOLS(f.v.order);
		f.v.poly.pp.fit(x,pp_f',w); 
	case {'direct'} % direct prediction
		f.v.meas_f = f.v.meas;
		if (f.nf.pre(1)>0)
			win     = kaiserwin(1:f.nf.pre(1));
			f.v.meas_f = wmedfilt(win,f.v.meas_f,1);
		end
		if (f.nf.pre(2)>0)
			win     = kaiserwin(1:f.nf.pre(2));
			f.v.meas_f = wmedfilt(win,f.v.meas_f',0)';
		end
		f.v.poly = PolyOLS(f.v.order);
		f.v.poly.fit(x,f.v.meas_f',w);
	otherwise
		error('here');
	end
	f.v.predict = @predict;
	f.v.pred    = f.v.predict(zs);
	% residual of predicted profile for reference profile
	f.v.res     = f.v.pred - f.v.meas;
	% residual by linearisation
	% f.v.res_lin = res_ln_z0.*df_dln_z0;

	% rmse 
	f.v.s      = sqrt(f.v.scale)*wrms(flat(f.v.weight(f.inner,:)),flat(f.v.res(f.inner,:)))
	% rmse locally along profile
	f.v.rmse = sqrt(f.v.scale)*wrms(f.v.weight,f.v.res,2);

	f.v           = f.profile_prediction_error(f.v,false,Q);

	if (f.nf.post>0)
		win     = kaiserwin(1:f.nf.post);
		f.v.rmse = sqrt(wmedfilt(win,f.v.rmse.^2,1));
	end	

	% TODO put into own file
	function pred = predict(x)
		switch (lower(f.v.mode)) % prediction from roughness lenght and wake parameter
		case {'profile'}
	                us_    = ones(size(us,1),length(x));                                           
			%us_ = us;
			ln_z0_ = f.v.poly.ln_z0.predict(x)';
			pp_    = f.v.poly.pp.predict(x)';
			zb_    = median(zb,2);
			
			% flow depth (h)
			h_ = bsxfun(@minus,rvec(x),zb_);
	
			% velocity at instrument depth (u)
			hi_ = zi - zb_;
			Si_ = bsxfun(@times,hi_,1./h_);
			%u  = l.predict(flat(Si)',flat(h)')';
			
			l = Log_profile_with_wake();
			l.shear_velocity(flat(us_));
			%l.roughness_length(flat(ln_z0_));
			l.ln_z0(flat(ln_z0_));
			l.perturbation_parameter(flat(pp_));

			% vertical profile at instrument level and partial derivative with respect to ln_z0
			% imaginary values for negative depths near the bank are cut off by inner selection
			% sensitivity
			%[f df_dln_z0] = log_profile(h_, 0, z_, ln_z0);
			%fp            = log_profile(h_, 0, z_, p_ln_z0');
	
		        u_    = l.u(flat(Si_)',flat(h_)')';                                    
	        	ubar_ = l.ubar(flat(h_).');

			u_    = real(reshape(u_,[],length(x)));
			ubar_ = real(reshape(ubar_,[],length(x)));
			pred  = (u_./ubar_);
		case {'direct'}
			pred = f.v.poly.predict(x)';
		otherwise
			error('here');
		end
	end % function pred

end % vertical_profile

